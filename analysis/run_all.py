#!/usr/bin/env python3
"""Reproducible local downstream analysis for the DrugRep HCC RNA-seq dataset."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "drugrep-matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom, spearmanr
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
PALETTE = {"NC": "#6B7280", "Sora": "#2563EB", "Spi": "#D69E2E", "Teg": "#C85C8E"}
QUANTIFIER_COLORS = {"RSEM": "#2563EB", "Salmon": "#D69E2E", "STAR": "#C85C8E"}


def now_iso() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def safe_tar_member(name: str) -> bool:
    member = Path(name)
    return not member.is_absolute() and ".." not in member.parts


def verify_internal_checksums(dataset_dir: Path) -> dict[str, Any]:
    sums_path = dataset_dir / "SHA256SUMS.txt"
    if not sums_path.exists():
        raise RuntimeError(f"Missing internal checksum file: {sums_path}")
    checked = 0
    for raw_line in sums_path.read_text(encoding="utf-8").splitlines():
        if not raw_line.strip():
            continue
        expected, rel_path = raw_line.split(maxsplit=1)
        rel_path = rel_path.strip().removeprefix("./")
        candidate = dataset_dir / rel_path
        if not candidate.is_file():
            raise RuntimeError(f"Checksum target is missing: {candidate}")
        actual = sha256_file(candidate)
        if actual != expected:
            raise RuntimeError(f"Checksum mismatch for {candidate}: {actual} != {expected}")
        checked += 1

    manifest_path = dataset_dir / "file_manifest.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t", header=None, names=["bytes", "path"])
    for row in manifest.itertuples(index=False):
        candidate = dataset_dir / str(row.path)
        if not candidate.is_file() or candidate.stat().st_size != int(row.bytes):
            raise RuntimeError(f"File-manifest mismatch: {candidate}")
    return {"checksums_verified": checked, "manifest_files_verified": int(len(manifest))}


def extract_archive(archive: Path, destination: Path) -> tuple[Path, dict[str, Any]]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        verification = verify_internal_checksums(destination)
        verification["reused_existing_extraction"] = True
        return destination, verification

    temp_root = Path(tempfile.mkdtemp(prefix=".extract-", dir=destination.parent))
    try:
        with tarfile.open(archive, "r:gz") as tf:
            members = tf.getmembers()
            if not members:
                raise RuntimeError("Archive is empty")
            for member in members:
                if not safe_tar_member(member.name):
                    raise RuntimeError(f"Unsafe archive member: {member.name}")
                if member.issym() or member.islnk():
                    raise RuntimeError(f"Links are not accepted in the archive: {member.name}")
                target = (temp_root / member.name).resolve()
                if temp_root.resolve() not in target.parents and target != temp_root.resolve():
                    raise RuntimeError(f"Archive member escapes extraction root: {member.name}")
            tf.extractall(temp_root, members=members)

        top_levels = sorted({Path(member.name).parts[0] for member in members if Path(member.name).parts})
        if len(top_levels) != 1:
            raise RuntimeError(f"Expected one archive root directory, found: {top_levels}")
        extracted = temp_root / top_levels[0]
        verification = verify_internal_checksums(extracted)
        extracted.rename(destination)
        verification["reused_existing_extraction"] = False
        return destination, verification
    finally:
        if temp_root.exists():
            shutil.rmtree(temp_root)


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def validate_inputs(dataset_dir: Path, config: dict[str, Any]) -> dict[str, Any]:
    groups = pd.read_csv(dataset_dir / "00_metadata_and_logs/group.list", sep="\t")
    expected_groups = config["groups"]
    if list(groups.columns) != ["sample", "group"]:
        raise RuntimeError("group.list must contain sample and group columns")
    if groups["sample"].duplicated().any():
        raise RuntimeError("Duplicate sample names in group.list")
    observed_counts = groups["group"].value_counts().to_dict()
    if set(observed_counts) != set(expected_groups) or any(observed_counts[g] != 3 for g in expected_groups):
        raise RuntimeError(f"Expected four groups with three replicates each, found {observed_counts}")

    expected_samples = groups["sample"].tolist()
    quantifier_profile: dict[str, Any] = {}
    raw_tables: dict[str, pd.DataFrame] = {}
    for quantifier in config["quantifiers"]:
        raw_path = dataset_dir / f"01_raw_counts/{quantifier}_gene_counts_matrix.txt"
        raw = read_table(raw_path)
        raw_tables[quantifier] = raw
        if raw.columns.tolist() != expected_samples:
            raise RuntimeError(f"Sample mismatch in {raw_path}")
        values = raw.to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            raise RuntimeError(f"Counts must be finite and nonnegative in {raw_path}")
        quantifier_profile[quantifier] = {
            "genes": int(raw.shape[0]),
            "samples": int(raw.shape[1]),
            "all_zero_genes": int((raw.sum(axis=1) == 0).sum()),
            "library_sizes": {col: int(raw[col].sum()) for col in raw.columns},
        }

        norm_path = dataset_dir / f"02_DESeq2_normalized/{quantifier}_DESeq2_normalize_matrix.txt"
        norm = read_table(norm_path)
        if norm.columns.tolist() != expected_samples or (norm.to_numpy() < 0).any():
            raise RuntimeError(f"Invalid normalized matrix: {norm_path}")

    deseq_columns = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
    for quantifier in config["quantifiers"]:
        for contrast in config["contrasts"]:
            de_path = dataset_dir / f"03_DESeq2_full/{quantifier}_DESeq2_{contrast}.txt"
            de = read_table(de_path)
            if not deseq_columns.issubset(de.columns):
                raise RuntimeError(f"DESeq2 schema mismatch: {de_path}")

    return {
        "samples": expected_samples,
        "group_counts": {key: int(value) for key, value in observed_counts.items()},
        "quantifier_profile": quantifier_profile,
    }


def prepare_idep_files(dataset_dir: Path, output_dir: Path, config: dict[str, Any]) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    quantifier = config["primary_quantifier"]
    raw = read_table(dataset_dir / f"01_raw_counts/{quantifier}_gene_counts_matrix.txt")
    raw.index.name = "gene_id"
    counts_path = output_dir / f"{quantifier}_raw_counts_idep.tsv"
    raw.to_csv(counts_path, sep="\t")
    groups = pd.read_csv(dataset_dir / "00_metadata_and_logs/group.list", sep="\t")
    design = groups.rename(columns={"sample": "Sample", "group": "Treatment"})
    design_path = output_dir / "sample_design_idep.tsv"
    design.to_csv(design_path, sep="\t", index=False)
    return [counts_path, design_path]


def write_preanalysis_inventory(dataset_dir: Path, table_dir: Path) -> Path:
    rows = [
        {
            "stage": "FASTQ quality control",
            "mentor_log_status": "Completed",
            "files_in_core_archive": "Pipeline/cluster logs only",
            "assessment": "Raw FastQC/MultiQC artifacts were not supplied; cannot independently audit read-level QC.",
        },
        {
            "stage": "Alignment",
            "mentor_log_status": "Completed",
            "files_in_core_archive": "No BAM or alignment summaries",
            "assessment": "Alignment was logged as complete but mapping metrics are unavailable in this archive.",
        },
        {
            "stage": "Expression quantification",
            "mentor_log_status": "Completed",
            "files_in_core_archive": "RSEM, Salmon, and STAR raw and normalized matrices",
            "assessment": "Core matrices are present and pass dimensions, sample, nonnegativity, manifest, and checksum checks.",
        },
        {
            "stage": "Differential expression",
            "mentor_log_status": "Completed",
            "files_in_core_archive": "9 DESeq2 full, 9 filtered, and 18 edgeR/limma tables",
            "assessment": "All archived DEG sets are exactly reproduced from the full tables at |log2FC| >= 1 and FDR < 0.05.",
        },
        {
            "stage": "Gene annotation post-processing",
            "mentor_log_status": "Failed: unexpected else",
            "files_in_core_archive": "Combined annotation CSVs intentionally excluded",
            "assessment": "Core expression/DE results remain usable; the local workflow performs an independent versioned gene-ID mapping.",
        },
    ]
    path = table_dir / "preanalysis_inventory.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def biomart_xml() -> str:
    return """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<!DOCTYPE Query>
<Query virtualSchemaName=\"default\" formatter=\"TSV\" header=\"0\" uniqueRows=\"1\" datasetConfigVersion=\"0.6\">
  <Dataset name=\"hsapiens_gene_ensembl\" interface=\"default\">
    <Attribute name=\"ensembl_gene_id\"/>
    <Attribute name=\"external_gene_name\"/>
    <Attribute name=\"entrezgene_id\"/>
    <Attribute name=\"gene_biotype\"/>
  </Dataset>
</Query>"""


def download_gene_map(path: Path) -> tuple[str, pd.DataFrame]:
    endpoints = [
        "https://www.ensembl.org/biomart/martservice",
        "https://useast.ensembl.org/biomart/martservice",
    ]
    last_error: Exception | None = None
    for endpoint in endpoints:
        for method in ("GET", "POST"):
            try:
                if method == "GET":
                    url = endpoint + "?" + urllib.parse.urlencode({"query": biomart_xml()})
                    request = urllib.request.Request(url, headers={"User-Agent": "DrugRep/1.0"})
                else:
                    body = urllib.parse.urlencode({"query": biomart_xml()}).encode("utf-8")
                    request = urllib.request.Request(endpoint, data=body, headers={"User-Agent": "DrugRep/1.0"})
                with urllib.request.urlopen(request, timeout=240) as response:
                    content = response.read()
                if content.startswith(b"Query ERROR") or len(content) < 1000:
                    raise RuntimeError(content[:500].decode("utf-8", errors="replace"))
                path.write_bytes(content)
                frame = pd.read_csv(
                    path,
                    sep="\t",
                    header=None,
                    names=["ensembl_gene_id", "gene_symbol", "entrez_gene_id", "gene_biotype"],
                    dtype=str,
                )
                frame = frame.dropna(subset=["ensembl_gene_id"]).drop_duplicates()
                frame.to_csv(path, sep="\t", index=False)
                return f"{endpoint} ({method})", frame
            except Exception as exc:  # try another method or the documented mirror
                last_error = exc
    # BioMart occasionally blocks automated requests. Fall back to a transient
    # official GENCODE GTF download and retain only the compact gene map.
    gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(prefix="gencode-", suffix=".gtf.gz", delete=False) as temp_handle:
            temp_path = Path(temp_handle.name)
            request = urllib.request.Request(gencode_url, headers={"User-Agent": "DrugRep/1.0"})
            with urllib.request.urlopen(request, timeout=600) as response:
                shutil.copyfileobj(response, temp_handle)
        attribute_pattern = re.compile(r'(gene_id|gene_name|gene_type) "([^"]+)"')
        rows: list[dict[str, str]] = []
        with gzip.open(temp_path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "gene":
                    continue
                attributes = dict(attribute_pattern.findall(fields[8]))
                gene_id = attributes.get("gene_id", "").split(".", 1)[0]
                if gene_id:
                    rows.append(
                        {
                            "ensembl_gene_id": gene_id,
                            "gene_symbol": attributes.get("gene_name", ""),
                            "entrez_gene_id": "",
                            "gene_biotype": attributes.get("gene_type", ""),
                        }
                    )
        frame = pd.DataFrame(rows).drop_duplicates("ensembl_gene_id", keep="first")
        if len(frame) < 50000:
            raise RuntimeError(f"GENCODE mapping unexpectedly contained only {len(frame)} genes")
        frame.to_csv(path, sep="\t", index=False)
        return "GENCODE release 49 (Ensembl 115) gene annotation", frame
    except Exception as exc:
        raise RuntimeError(f"Unable to download Ensembl BioMart or GENCODE mapping: BioMart={last_error}; GENCODE={exc}") from exc
    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()


def write_gmt(path: Path, gene_sets: dict[str, list[str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for term in sorted(gene_sets):
            genes = sorted({str(gene).strip().upper() for gene in gene_sets[term] if str(gene).strip()})
            if genes:
                handle.write("\t".join([term, "Enrichr"] + genes) + "\n")


def read_gmt(path: Path) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                result[fields[0]] = [gene.upper() for gene in fields[2:] if gene]
    return result


def library_prefix(name: str) -> str:
    for prefix in ("MSigDB_Hallmark", "GO_Biological_Process", "KEGG", "Reactome"):
        if name.startswith(prefix):
            return prefix
    return name


def library_sort_key(name: str) -> tuple[int, str]:
    years = [int(value) for value in re.findall(r"(?:19|20)\d{2}", name)]
    return (max(years) if years else 0, name)


def prepare_references(
    reference_dir: Path,
    config: dict[str, Any],
    offline: bool,
    refresh: bool,
    warnings: list[str],
) -> tuple[pd.DataFrame, dict[str, Path], dict[str, Any]]:
    reference_dir.mkdir(parents=True, exist_ok=True)
    gmt_dir = reference_dir / "gmt"
    gmt_dir.mkdir(parents=True, exist_ok=True)
    mapping_path = reference_dir / "ensembl_human_gene_map.tsv"
    manifest_path = reference_dir / "reference_manifest.json"

    if offline and refresh:
        raise RuntimeError("--offline and --refresh-references cannot be combined")
    if refresh and not offline:
        for candidate in [mapping_path, manifest_path]:
            if candidate.exists():
                candidate.unlink()
        for candidate in gmt_dir.glob("*.gmt"):
            candidate.unlink()

    if mapping_path.exists():
        mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
        mapping_source = "cached Ensembl BioMart snapshot"
    elif offline:
        raise RuntimeError(f"Offline mode requires {mapping_path}")
    else:
        mapping_source, mapping = download_gene_map(mapping_path)

    try:
        import gseapy as gp
    except ImportError as exc:
        raise RuntimeError("GSEApy is required; install requirements.txt") from exc

    selected_paths: dict[str, Path] = {}
    selected_names: dict[str, str] = {}
    available: list[str] | None = None
    for requested in config["gene_sets"]:
        exact_path = gmt_dir / f"{requested}.gmt"
        if exact_path.exists():
            selected_paths[requested] = exact_path
            selected_names[requested] = requested
            continue
        if offline:
            raise RuntimeError(f"Offline mode requires cached gene set: {exact_path}")
        if available is None:
            available = gp.get_library_name(organism="Human")
        actual = requested
        if requested not in available:
            prefix = library_prefix(requested)
            candidates = sorted((name for name in available if name.startswith(prefix)), key=library_sort_key)
            if not candidates:
                raise RuntimeError(f"No Enrichr library matches {requested}")
            actual = candidates[-1]
            warnings.append(f"Requested {requested}; cached latest matching Enrichr library {actual}.")
        gene_sets = gp.get_library(name=actual, organism="Human")
        write_gmt(exact_path, gene_sets)
        selected_paths[requested] = exact_path
        selected_names[requested] = actual

    references: list[dict[str, Any]] = [
        {
            "type": "gene_mapping",
            "name": "Ensembl human gene mapping",
            "source": mapping_source,
            "retrieved_or_verified_at": now_iso(),
            "path": relative(mapping_path),
            "rows": int(len(mapping)),
            "sha256": sha256_file(mapping_path),
            "license_note": "GENCODE/Ensembl annotation; review upstream data-use and attribution terms before redistribution.",
        }
    ]
    for requested, path in selected_paths.items():
        references.append(
            {
                "type": "gene_set_library",
                "requested_name": requested,
                "selected_name": selected_names[requested],
                "source": "Enrichr library service via GSEApy",
                "retrieved_or_verified_at": now_iso(),
                "path": relative(path),
                "gene_set_count": len(read_gmt(path)),
                "sha256": sha256_file(path),
                "license_note": "Terms inherit from each upstream gene-set source; review before redistribution.",
            }
        )
    reference_manifest = {"created_at": now_iso(), "references": references}
    manifest_path.write_text(json.dumps(reference_manifest, indent=2), encoding="utf-8")
    return mapping, selected_paths, reference_manifest


def collapse_mapping(mapping: pd.DataFrame) -> pd.DataFrame:
    required = {"ensembl_gene_id", "gene_symbol"}
    if not required.issubset(mapping.columns):
        raise RuntimeError(f"Mapping must contain {required}")
    clean = mapping.copy()
    clean["gene_symbol"] = clean["gene_symbol"].fillna("").astype(str).str.strip().str.upper()
    clean = clean[clean["gene_symbol"] != ""].sort_values(["ensembl_gene_id", "gene_symbol"])
    return clean.drop_duplicates("ensembl_gene_id", keep="first")


def build_ranked_table(de: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    clean_map = collapse_mapping(mapping).set_index("ensembl_gene_id")
    ranked = de[["stat", "log2FoldChange", "padj"]].copy()
    ranked.index = ranked.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    ranked["ensembl_gene_id"] = ranked.index
    ranked = ranked.reset_index(drop=True)
    ranked["gene_symbol"] = ranked["ensembl_gene_id"].map(clean_map["gene_symbol"])
    ranked = ranked.dropna(subset=["stat", "gene_symbol"]).copy()
    ranked["abs_stat"] = ranked["stat"].abs()
    ranked = ranked.sort_values(["gene_symbol", "abs_stat", "ensembl_gene_id"], ascending=[True, False, True])
    ranked = ranked.drop_duplicates("gene_symbol", keep="first")
    ranked = ranked.sort_values(["stat", "gene_symbol"], ascending=[False, True]).reset_index(drop=True)
    # Break exact statistic ties deterministically by moving tied followers one
    # representable float toward -infinity. This preserves practical scale,
    # direction, and the secondary gene-symbol order established above.
    scores = ranked["stat"].to_numpy(dtype=float).copy()
    for index in range(1, len(scores)):
        if scores[index] >= scores[index - 1]:
            scores[index] = np.nextafter(scores[index - 1], -np.inf)
    ranked["rank_score"] = scores
    return ranked[["ensembl_gene_id", "gene_symbol", "rank_score", "stat", "log2FoldChange", "padj"]]


def deg_mask(frame: pd.DataFrame, config: dict[str, Any], lfc_col: str, fdr_col: str) -> pd.Series:
    return frame[fdr_col].notna() & (frame[fdr_col] < config["deg"]["padj"]) & (
        frame[lfc_col].abs() >= config["deg"]["abs_log2fc"]
    )


def reproduce_deg_tables(dataset_dir: Path, config: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for quantifier in config["quantifiers"]:
        for contrast in config["contrasts"]:
            full = read_table(dataset_dir / f"03_DESeq2_full/{quantifier}_DESeq2_{contrast}.txt")
            archived = read_table(
                dataset_dir / f"04_DESeq2_significant/{quantifier}_DESeq2_{contrast}_filter_1_0.05.txt"
            )
            selected = full.loc[deg_mask(full, config, "log2FoldChange", "padj")]
            if set(selected.index) != set(archived.index):
                missing = sorted(set(archived.index) - set(selected.index))[:5]
                extra = sorted(set(selected.index) - set(archived.index))[:5]
                raise RuntimeError(
                    f"DEG regeneration mismatch for {quantifier} {contrast}; missing={missing}, extra={extra}"
                )
            rows.append(
                {
                    "quantifier": quantifier,
                    "contrast": contrast,
                    "tested_genes": int(full["stat"].notna().sum()),
                    "padj_missing": int(full["padj"].isna().sum()),
                    "deg_total": int(len(selected)),
                    "deg_up": int((selected["log2FoldChange"] > 0).sum()),
                    "deg_down": int((selected["log2FoldChange"] < 0).sum()),
                    "archived_deg_reproduced": True,
                }
            )
    return pd.DataFrame(rows)


def spearman_value(x: pd.Series, y: pd.Series) -> float:
    valid = x.notna() & y.notna()
    if valid.sum() < 3:
        return math.nan
    return float(spearmanr(x[valid], y[valid]).statistic)


def jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    a, b = set(left), set(right)
    return len(a & b) / len(a | b) if a or b else math.nan


def analyze_de(
    dataset_dir: Path,
    mapping: pd.DataFrame,
    processed_dir: Path,
    table_dir: Path,
    config: dict[str, Any],
) -> dict[str, Any]:
    de_dir = processed_dir / "differential_expression"
    rank_dir = processed_dir / "ranked"
    de_dir.mkdir(parents=True, exist_ok=True)
    rank_dir.mkdir(parents=True, exist_ok=True)
    clean_map = collapse_mapping(mapping).set_index("ensembl_gene_id")

    all_de: dict[tuple[str, str], pd.DataFrame] = {}
    ranked_tables: dict[tuple[str, str], pd.DataFrame] = {}
    mapping_rows: list[dict[str, Any]] = []
    top_rows: list[pd.DataFrame] = []
    for quantifier in config["quantifiers"]:
        for contrast in config["contrasts"]:
            frame = read_table(dataset_dir / f"03_DESeq2_full/{quantifier}_DESeq2_{contrast}.txt")
            frame.index = frame.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
            frame.index.name = "ensembl_gene_id"
            frame["gene_symbol"] = frame.index.map(clean_map["gene_symbol"])
            all_de[(quantifier, contrast)] = frame
            frame.to_csv(de_dir / f"{quantifier}_DESeq2_{contrast}_mapped.tsv", sep="\t")
            ranked = build_ranked_table(frame, mapping)
            ranked_tables[(quantifier, contrast)] = ranked
            ranked[["gene_symbol", "rank_score"]].to_csv(
                rank_dir / f"{quantifier}_{contrast}.rnk", sep="\t", index=False, header=False
            )
            tested_ids = frame.index[frame["stat"].notna()]
            mapped_tested = frame.loc[tested_ids, "gene_symbol"].notna().sum()
            mapping_rows.append(
                {
                    "quantifier": quantifier,
                    "contrast": contrast,
                    "tested_genes": int(len(tested_ids)),
                    "mapped_tested_genes": int(mapped_tested),
                    "mapping_coverage": float(mapped_tested / len(tested_ids)),
                    "unique_ranked_symbols": int(len(ranked)),
                }
            )
            if quantifier == config["primary_quantifier"]:
                ranked_top = ranked.reindex(ranked["stat"].abs().sort_values(ascending=False).index).head(20).copy()
                ranked_top.insert(0, "contrast", contrast)
                top_rows.append(ranked_top)

    mapping_summary = pd.DataFrame(mapping_rows)
    mapping_summary.to_csv(table_dir / "mapping_summary.tsv", sep="\t", index=False)
    if mapping_summary["mapping_coverage"].min() < config["mapping"]["minimum_coverage"]:
        raise RuntimeError(
            f"Gene-symbol mapping coverage {mapping_summary['mapping_coverage'].min():.1%} is below "
            f"{config['mapping']['minimum_coverage']:.1%}"
        )
    pd.concat(top_rows, ignore_index=True).to_csv(table_dir / "top_genes_by_wald_stat.tsv", sep="\t", index=False)

    robustness: list[dict[str, Any]] = []
    primary = config["primary_quantifier"]
    for contrast in config["contrasts"]:
        base = all_de[(primary, contrast)]
        base_deg = set(base.index[deg_mask(base, config, "log2FoldChange", "padj")])
        for quantifier in config["quantifiers"]:
            if quantifier == primary:
                continue
            other = all_de[(quantifier, contrast)]
            common = base.index.intersection(other.index)
            other_deg = set(other.index[deg_mask(other, config, "log2FoldChange", "padj")])
            both_sig = sorted(base_deg & other_deg)
            sign_agreement = (
                float((np.sign(base.loc[both_sig, "log2FoldChange"]) == np.sign(other.loc[both_sig, "log2FoldChange"])).mean())
                if both_sig
                else math.nan
            )
            robustness.append(
                {
                    "contrast": contrast,
                    "comparison": f"{primary} DESeq2 vs {quantifier} DESeq2",
                    "log2fc_spearman": spearman_value(
                        base.loc[common, "log2FoldChange"], other.loc[common, "log2FoldChange"]
                    ),
                    "stat_spearman": spearman_value(base.loc[common, "stat"], other.loc[common, "stat"]),
                    "deg_jaccard": jaccard(base_deg, other_deg),
                    "shared_deg_sign_agreement": sign_agreement,
                    "common_genes": int(len(common)),
                }
            )

        for method, lfc_col, fdr_col, stat_col in [
            ("edgeR", "logFC", "FDR", "F"),
            ("Limma", "logFC", "adj.P.Val", "t"),
        ]:
            other = read_table(dataset_dir / f"05_edgeR_Limma_sensitivity/{primary}_{method}_{contrast}.txt")
            common = base.index.intersection(other.index)
            other_deg = set(other.index[deg_mask(other, config, lfc_col, fdr_col)])
            both_sig = sorted(base_deg & other_deg)
            robustness.append(
                {
                    "contrast": contrast,
                    "comparison": f"{primary} DESeq2 vs {primary} {method}",
                    "log2fc_spearman": spearman_value(base.loc[common, "log2FoldChange"], other.loc[common, lfc_col]),
                    "stat_spearman": (
                        spearman_value(base.loc[common, "stat"], other.loc[common, stat_col])
                        if method == "Limma"
                        else math.nan
                    ),
                    "deg_jaccard": jaccard(base_deg, other_deg),
                    "shared_deg_sign_agreement": (
                        float((np.sign(base.loc[both_sig, "log2FoldChange"]) == np.sign(other.loc[both_sig, lfc_col])).mean())
                        if both_sig
                        else math.nan
                    ),
                    "common_genes": int(len(common)),
                }
            )
    robustness_frame = pd.DataFrame(robustness)
    robustness_frame.to_csv(table_dir / "robustness_summary.tsv", sep="\t", index=False)

    similarity: list[dict[str, Any]] = []
    contrasts = list(config["contrasts"])
    for i, left_name in enumerate(contrasts):
        for right_name in contrasts[i + 1 :]:
            left = all_de[(primary, left_name)]
            right = all_de[(primary, right_name)]
            common = left.index.intersection(right.index)
            left_deg = set(left.index[deg_mask(left, config, "log2FoldChange", "padj")])
            right_deg = set(right.index[deg_mask(right, config, "log2FoldChange", "padj")])
            shared = sorted(left_deg & right_deg)
            similarity.append(
                {
                    "contrast_a": left_name,
                    "contrast_b": right_name,
                    "stat_spearman": spearman_value(left.loc[common, "stat"], right.loc[common, "stat"]),
                    "log2fc_spearman": spearman_value(
                        left.loc[common, "log2FoldChange"], right.loc[common, "log2FoldChange"]
                    ),
                    "deg_jaccard": jaccard(left_deg, right_deg),
                    "shared_deg": len(shared),
                    "shared_deg_sign_agreement": (
                        float((np.sign(left.loc[shared, "log2FoldChange"]) == np.sign(right.loc[shared, "log2FoldChange"])).mean())
                        if shared
                        else math.nan
                    ),
                }
            )
    similarity_frame = pd.DataFrame(similarity)
    similarity_frame.to_csv(table_dir / "treatment_similarity.tsv", sep="\t", index=False)
    return {
        "all_de": all_de,
        "ranked": ranked_tables,
        "mapping_summary": mapping_summary,
        "robustness": robustness_frame,
        "similarity": similarity_frame,
    }


def run_qc(dataset_dir: Path, table_dir: Path, config: dict[str, Any]) -> dict[str, Any]:
    primary = config["primary_quantifier"]
    raw = read_table(dataset_dir / f"01_raw_counts/{primary}_gene_counts_matrix.txt")
    norm = read_table(dataset_dir / f"02_DESeq2_normalized/{primary}_DESeq2_normalize_matrix.txt")
    groups = pd.read_csv(dataset_dir / "00_metadata_and_logs/group.list", sep="\t").set_index("sample")
    log_expression = np.log2(norm + 1)
    library_sizes = raw.sum(axis=0).rename("library_size").to_frame()
    library_sizes["million_reads"] = library_sizes["library_size"] / 1e6
    library_sizes["group"] = groups.loc[library_sizes.index, "group"]
    library_sizes.index.name = "sample"
    library_sizes.to_csv(table_dir / "library_sizes.tsv", sep="\t")

    variances = log_expression.var(axis=1).sort_values(ascending=False)
    top = log_expression.loc[variances.head(config["qc"]["top_variable_genes"]).index]
    pca = PCA(n_components=min(4, top.shape[1] - 1))
    scores = pca.fit_transform(top.T)
    score_frame = pd.DataFrame(scores, index=top.columns, columns=[f"PC{i + 1}" for i in range(scores.shape[1])])
    score_frame["group"] = groups.loc[score_frame.index, "group"]
    score_frame.index.name = "sample"
    score_frame.to_csv(table_dir / "pca_scores.tsv", sep="\t")

    correlations = log_expression.corr(method="pearson")
    correlations.to_csv(table_dir / "sample_correlations.tsv", sep="\t")
    outliers: list[dict[str, Any]] = []
    floor = config["qc"]["outlier_correlation_floor"]
    for sample in correlations.index:
        group = groups.loc[sample, "group"]
        peers = [name for name in correlations.columns if groups.loc[name, "group"] == group and name != sample]
        mean_peer_correlation = float(correlations.loc[sample, peers].mean())
        outliers.append(
            {
                "sample": sample,
                "group": group,
                "mean_within_group_correlation": mean_peer_correlation,
                "flagged_below_floor": bool(mean_peer_correlation < floor),
            }
        )
    outlier_frame = pd.DataFrame(outliers)
    outlier_frame.to_csv(table_dir / "sample_outlier_flags.tsv", sep="\t", index=False)
    return {
        "raw": raw,
        "norm": norm,
        "log_expression": log_expression,
        "groups": groups,
        "library_sizes": library_sizes,
        "pca_scores": score_frame,
        "pca_variance": pca.explained_variance_ratio_,
        "correlations": correlations,
        "top_expression": top,
        "outliers": outlier_frame,
    }


def normalize_gsea_result(result: pd.DataFrame) -> pd.DataFrame:
    frame = result.copy()
    if "Term" not in frame.columns:
        frame = frame.reset_index(drop=False)
        first = frame.columns[0]
        frame = frame.rename(columns={first: "Term"})
    else:
        # GSEApy may retain a process-local integer index whose ordering is not
        # semantically meaningful. Never serialize it into reproducible output.
        frame = frame.reset_index(drop=True)
    aliases = {
        "Name": "gsea_name",
        "Term": "term",
        "ES": "es",
        "NES": "nes",
        "NOM p-val": "nominal_p",
        "FDR q-val": "fdr",
        "FWER p-val": "fwer",
        "Tag %": "tag_fraction",
        "Gene %": "gene_fraction",
        "Lead_genes": "leading_edge_genes",
    }
    frame = frame.rename(columns={key: value for key, value in aliases.items() if key in frame.columns})
    for column in ["es", "nes", "nominal_p", "fdr", "fwer"]:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame.sort_values(["term"], kind="mergesort").reset_index(drop=True)


def run_ora(
    gene_list: set[str],
    background: set[str],
    gene_sets: dict[str, list[str]],
    min_size: int,
    max_size: int,
) -> pd.DataFrame:
    background = {gene.upper() for gene in background}
    selected = {gene.upper() for gene in gene_list} & background
    rows: list[dict[str, Any]] = []
    for term, members_raw in gene_sets.items():
        members = set(members_raw) & background
        if not min_size <= len(members) <= max_size:
            continue
        overlap = selected & members
        if len(overlap) < 2:
            continue
        pvalue = float(hypergeom.sf(len(overlap) - 1, len(background), len(members), len(selected)))
        rows.append(
            {
                "term": term,
                "overlap": len(overlap),
                "query_size": len(selected),
                "gene_set_size": len(members),
                "background_size": len(background),
                "pvalue": pvalue,
                "overlap_genes": ";".join(sorted(overlap)),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(
            columns=["term", "overlap", "query_size", "gene_set_size", "background_size", "pvalue", "fdr", "overlap_genes"]
        )
    frame["fdr"] = multipletests(frame["pvalue"].to_numpy(), method="fdr_bh")[1]
    return frame.sort_values(["fdr", "pvalue", "term"]).reset_index(drop=True)


def run_enrichment(
    de_results: dict[str, Any],
    gmt_paths: dict[str, Path],
    table_dir: Path,
    config: dict[str, Any],
) -> dict[str, pd.DataFrame]:
    try:
        import gseapy as gp
    except ImportError as exc:
        raise RuntimeError("GSEApy is required; install requirements.txt") from exc

    gsea_frames: list[pd.DataFrame] = []
    ora_frames: list[pd.DataFrame] = []
    primary = config["primary_quantifier"]
    gsea_config = config["gsea"]
    hallmark = config["hallmark_library"]

    jobs: list[tuple[str, str, str]] = []
    for contrast in config["contrasts"]:
        for library in config["gene_sets"]:
            jobs.append((primary, contrast, library))
    for quantifier in config["quantifiers"]:
        if quantifier != primary:
            for contrast in config["contrasts"]:
                jobs.append((quantifier, contrast, hallmark))

    for quantifier, contrast, library in jobs:
        ranked = de_results["ranked"][(quantifier, contrast)]
        rnk = ranked[["gene_symbol", "rank_score"]].copy()
        result = gp.prerank(
            rnk=rnk,
            gene_sets=str(gmt_paths[library]),
            threads=gsea_config["threads"],
            min_size=gsea_config["min_size"],
            max_size=gsea_config["max_size"],
            permutation_num=gsea_config["permutations"],
            outdir=None,
            seed=gsea_config["seed"],
            verbose=False,
        )
        frame = normalize_gsea_result(result.res2d)
        frame.insert(0, "library", library)
        frame.insert(0, "contrast", contrast)
        frame.insert(0, "quantifier", quantifier)
        gsea_frames.append(frame)

    gsea_all = pd.concat(gsea_frames, ignore_index=True)
    gsea_all = gsea_all.sort_values(
        ["quantifier", "contrast", "library", "term"], kind="mergesort"
    ).reset_index(drop=True)
    gsea_all.to_csv(
        table_dir / "gsea_all.tsv", sep="\t", index=False, float_format="%.12g"
    )

    for contrast in config["contrasts"]:
        frame = de_results["all_de"][(primary, contrast)]
        mapped = de_results["ranked"][(primary, contrast)]
        background = set(mapped["gene_symbol"])
        significant = frame.loc[deg_mask(frame, config, "log2FoldChange", "padj")].copy()
        significant["gene_symbol"] = significant["gene_symbol"].astype("string").str.upper()
        for direction, subset in [
            ("up", significant[significant["log2FoldChange"] > 0]),
            ("down", significant[significant["log2FoldChange"] < 0]),
        ]:
            gene_list = set(subset["gene_symbol"].dropna())
            for library, path in gmt_paths.items():
                ora = run_ora(
                    gene_list,
                    background,
                    read_gmt(path),
                    gsea_config["min_size"],
                    gsea_config["max_size"],
                )
                ora.insert(0, "library", library)
                ora.insert(0, "direction", direction)
                ora.insert(0, "contrast", contrast)
                ora_frames.append(ora)
    ora_all = pd.concat(ora_frames, ignore_index=True)
    ora_all.to_csv(table_dir / "ora_all.tsv", sep="\t", index=False)

    sensitivity_rows: list[dict[str, Any]] = []
    for contrast in config["contrasts"]:
        base = gsea_all[
            (gsea_all["quantifier"] == primary)
            & (gsea_all["contrast"] == contrast)
            & (gsea_all["library"] == hallmark)
        ].set_index("term")
        for quantifier in config["quantifiers"]:
            if quantifier == primary:
                continue
            other = gsea_all[
                (gsea_all["quantifier"] == quantifier)
                & (gsea_all["contrast"] == contrast)
                & (gsea_all["library"] == hallmark)
            ].set_index("term")
            common = base.index.intersection(other.index)
            exploratory = common[
                (base.loc[common, "fdr"] < gsea_config["exploratory_fdr"])
                | (other.loc[common, "fdr"] < gsea_config["exploratory_fdr"])
            ]
            sensitivity_rows.append(
                {
                    "contrast": contrast,
                    "comparison": f"{primary} vs {quantifier}",
                    "nes_spearman": spearman_value(base.loc[common, "nes"], other.loc[common, "nes"]),
                    "exploratory_pathways": int(len(exploratory)),
                    "exploratory_sign_agreement": (
                        float((np.sign(base.loc[exploratory, "nes"]) == np.sign(other.loc[exploratory, "nes"])).mean())
                        if len(exploratory)
                        else math.nan
                    ),
                }
            )
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(table_dir / "hallmark_sensitivity.tsv", sep="\t", index=False)

    pathway_similarity: list[dict[str, Any]] = []
    primary_gsea = gsea_all[gsea_all["quantifier"] == primary]
    contrasts = list(config["contrasts"])
    for i, left in enumerate(contrasts):
        for right in contrasts[i + 1 :]:
            a = primary_gsea[primary_gsea["contrast"] == left].set_index(["library", "term"])
            b = primary_gsea[primary_gsea["contrast"] == right].set_index(["library", "term"])
            common = a.index.intersection(b.index)
            pathway_similarity.append(
                {
                    "contrast_a": left,
                    "contrast_b": right,
                    "nes_spearman": spearman_value(a.loc[common, "nes"], b.loc[common, "nes"]),
                    "common_pathways": int(len(common)),
                }
            )
    pathway_similarity_frame = pd.DataFrame(pathway_similarity)
    pathway_similarity_frame.to_csv(table_dir / "pathway_similarity.tsv", sep="\t", index=False)

    leading = gsea_all[
        (gsea_all["quantifier"] == primary) & (gsea_all["fdr"] < gsea_config["main_fdr"])
    ].copy()
    columns = [
        "contrast",
        "library",
        "term",
        "nes",
        "nominal_p",
        "fdr",
        "leading_edge_genes",
    ]
    leading[[column for column in columns if column in leading.columns]].to_csv(
        table_dir / "gsea_significant_leading_edges.tsv",
        sep="\t",
        index=False,
        float_format="%.12g",
    )
    return {
        "gsea": gsea_all,
        "ora": ora_all,
        "hallmark_sensitivity": sensitivity,
        "pathway_similarity": pathway_similarity_frame,
    }


def save_figure(fig: plt.Figure, base_path: Path) -> None:
    fig.savefig(base_path.with_suffix(".png"), dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(base_path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def group_from_sample(sample: str) -> str:
    return sample.split(".", 1)[0]


def plot_qc(qc: dict[str, Any], figure_dir: Path) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    library = qc["library_sizes"].reset_index()
    axes[0, 0].bar(
        library["sample"],
        library["million_reads"],
        color=[PALETTE[group] for group in library["group"]],
        edgecolor="#374151",
        linewidth=0.5,
    )
    axes[0, 0].set_title("RSEM library sizes")
    axes[0, 0].set_ylabel("Total assigned reads (millions)")
    axes[0, 0].tick_params(axis="x", rotation=45)

    scores = qc["pca_scores"].reset_index()
    for group, subset in scores.groupby("group", sort=False):
        axes[0, 1].scatter(subset["PC1"], subset["PC2"], s=70, label=group, color=PALETTE[group], edgecolor="white")
        for row in subset.itertuples(index=False):
            axes[0, 1].annotate(row.sample, (row.PC1, row.PC2), xytext=(4, 4), textcoords="offset points", fontsize=8)
    axes[0, 1].set_title("PCA of 500 most-variable genes")
    axes[0, 1].set_xlabel(f"PC1 ({qc['pca_variance'][0] * 100:.1f}% variance)")
    axes[0, 1].set_ylabel(f"PC2 ({qc['pca_variance'][1] * 100:.1f}% variance)")
    axes[0, 1].legend(frameon=False, ncol=2)

    sns.heatmap(qc["correlations"], cmap="Blues", vmin=0.96, vmax=1.0, square=True, ax=axes[1, 0], cbar_kws={"label": "Pearson r"})
    axes[1, 0].set_title("Sample correlation on log2 normalized counts")
    axes[1, 0].tick_params(axis="x", rotation=45)
    axes[1, 0].tick_params(axis="y", rotation=0)

    top50 = qc["top_expression"].head(50)
    z = top50.sub(top50.mean(axis=1), axis=0).div(top50.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    sns.heatmap(z, cmap="vlag", center=0, xticklabels=True, yticklabels=False, ax=axes[1, 1], cbar_kws={"label": "Gene-wise z-score"})
    axes[1, 1].set_title("Top 50 variable expression profiles")
    axes[1, 1].tick_params(axis="x", rotation=45)
    fig.suptitle("RNA-seq sample quality and expression structure", fontsize=15, fontweight="bold")
    save_figure(fig, figure_dir / "01_qc_overview")


def plot_deg(
    de_results: dict[str, Any],
    deg_summary: pd.DataFrame,
    figure_dir: Path,
    config: dict[str, Any],
) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.9)
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
    primary = config["primary_quantifier"]
    summary = deg_summary[deg_summary["quantifier"] == primary].copy()
    x = np.arange(len(summary))
    axes[0, 0].bar(x, summary["deg_up"], color="#D69E2E", label="Up")
    axes[0, 0].bar(x, summary["deg_down"], bottom=summary["deg_up"], color="#2563EB", label="Down")
    axes[0, 0].set_xticks(x, [name.replace("vsNC", " vs NC") for name in summary["contrast"]])
    axes[0, 0].set_ylabel("DEGs")
    axes[0, 0].set_title("RSEM DESeq2 DEG counts")
    axes[0, 0].legend(frameon=False)

    for ax, contrast in zip(axes.flat[1:], config["contrasts"]):
        frame = de_results["all_de"][(primary, contrast)].copy()
        frame["minus_log10_padj"] = -np.log10(frame["padj"].clip(lower=np.finfo(float).tiny))
        significant = deg_mask(frame, config, "log2FoldChange", "padj")
        colors = np.where(significant & (frame["log2FoldChange"] > 0), "#D69E2E", np.where(significant, "#2563EB", "#C7CDD4"))
        ax.scatter(frame["log2FoldChange"], frame["minus_log10_padj"], c=colors, s=8, alpha=0.65, linewidths=0)
        ax.axvline(config["deg"]["abs_log2fc"], color="#6B7280", linestyle="--", linewidth=0.8)
        ax.axvline(-config["deg"]["abs_log2fc"], color="#6B7280", linestyle="--", linewidth=0.8)
        ax.axhline(-math.log10(config["deg"]["padj"]), color="#6B7280", linestyle="--", linewidth=0.8)
        labels = frame.loc[significant].sort_values("padj").head(5)
        for gene_id, row in labels.iterrows():
            label = row["gene_symbol"] if pd.notna(row["gene_symbol"]) else gene_id
            ax.annotate(str(label), (row["log2FoldChange"], row["minus_log10_padj"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
        ax.set_title(contrast.replace("vsNC", " vs NC"))
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 adjusted p-value")
    fig.suptitle("Differential-expression response by treatment", fontsize=15, fontweight="bold")
    save_figure(fig, figure_dir / "02_differential_expression")


def plot_concordance(de_results: dict[str, Any], figure_dir: Path, config: dict[str, Any]) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.85)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    primary = config["primary_quantifier"]
    contrasts = list(config["contrasts"])
    labels = [name.replace("vsNC", "") for name in contrasts]
    stat_matrix = pd.DataFrame(np.eye(len(contrasts)), index=labels, columns=labels)
    jaccard_matrix = pd.DataFrame(np.eye(len(contrasts)), index=labels, columns=labels)
    for row in de_results["similarity"].itertuples(index=False):
        a, b = row.contrast_a.replace("vsNC", ""), row.contrast_b.replace("vsNC", "")
        stat_matrix.loc[a, b] = stat_matrix.loc[b, a] = row.stat_spearman
        jaccard_matrix.loc[a, b] = jaccard_matrix.loc[b, a] = row.deg_jaccard
    sns.heatmap(stat_matrix, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1, square=True, ax=axes[0, 0])
    axes[0, 0].set_title("Treatment similarity: Wald statistics")
    sns.heatmap(jaccard_matrix, annot=True, fmt=".2f", cmap="YlOrBr", vmin=0, vmax=1, square=True, ax=axes[0, 1])
    axes[0, 1].set_title("Treatment similarity: DEG Jaccard index")

    robustness = de_results["robustness"].copy()
    pivot = robustness.pivot(index="comparison", columns="contrast", values="log2fc_spearman")
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="Blues", vmin=0.7, vmax=1.0, ax=axes[1, 0], cbar_kws={"label": "Spearman ρ"})
    axes[1, 0].set_title("Cross-quantifier/method log2FC concordance")
    axes[1, 0].set_xlabel("")
    axes[1, 0].set_ylabel("")

    sora = de_results["all_de"][(primary, "SoravsNC")]
    for contrast, marker, color in [("SpivsNC", "o", PALETTE["Spi"]), ("TegvsNC", "^", PALETTE["Teg"])]:
        other = de_results["all_de"][(primary, contrast)]
        common = sora.index.intersection(other.index)
        axes[1, 1].scatter(
            sora.loc[common, "stat"],
            other.loc[common, "stat"],
            s=6,
            alpha=0.18,
            color=color,
            marker=marker,
            label=contrast.replace("vsNC", ""),
            linewidths=0,
        )
    axes[1, 1].axhline(0, color="#9CA3AF", linewidth=0.7)
    axes[1, 1].axvline(0, color="#9CA3AF", linewidth=0.7)
    axes[1, 1].set_xlabel("Sorafenib vs NC Wald statistic")
    axes[1, 1].set_ylabel("Candidate vs NC Wald statistic")
    axes[1, 1].set_title("Candidate perturbations versus sorafenib")
    axes[1, 1].legend(frameon=False)
    fig.suptitle("Treatment and analytical robustness", fontsize=15, fontweight="bold")
    save_figure(fig, figure_dir / "03_concordance")


def shorten_term(term: str, limit: int = 52) -> str:
    clean = re.sub(r"^(HALLMARK_|KEGG_|REACTOME_|GO[_ :]*)", "", str(term), flags=re.IGNORECASE)
    clean = clean.replace("_", " ").title()
    return clean if len(clean) <= limit else clean[: limit - 1] + "…"


def plot_pathways(enrichment: dict[str, pd.DataFrame], figure_dir: Path, config: dict[str, Any]) -> None:
    sns.set_theme(style="whitegrid", font_scale=0.8)
    fig, axes = plt.subplots(2, 2, figsize=(15, 12), constrained_layout=True)
    primary = config["primary_quantifier"]
    gsea = enrichment["gsea"]
    primary_gsea = gsea[gsea["quantifier"] == primary].copy()
    hallmark_gsea = primary_gsea[primary_gsea["library"] == config["hallmark_library"]].copy()
    hallmark_gsea["display_term"] = hallmark_gsea["term"].map(shorten_term)
    significant_terms = hallmark_gsea.loc[
        hallmark_gsea["fdr"] < config["gsea"]["main_fdr"], "display_term"
    ].unique()
    heat = hallmark_gsea[hallmark_gsea["display_term"].isin(significant_terms)].pivot_table(
        index="display_term", columns="contrast", values="nes", aggfunc="first"
    )
    heat = heat.reindex(heat.abs().max(axis=1).sort_values(ascending=False).head(20).index)
    sns.heatmap(heat, cmap="vlag", center=0, ax=axes[0, 0], cbar_kws={"label": "Normalized enrichment score"})
    axes[0, 0].set_title("Strongest significant Hallmark GSEA responses")
    axes[0, 0].set_xlabel("")
    axes[0, 0].set_ylabel("")

    for ax, contrast in zip([axes[0, 1], axes[1, 0]], ["SpivsNC", "TegvsNC"]):
        subset = hallmark_gsea[
            (hallmark_gsea["contrast"] == contrast)
            & (hallmark_gsea["fdr"] < config["gsea"]["exploratory_fdr"])
        ].copy()
        subset = subset.reindex(subset["nes"].abs().sort_values(ascending=False).head(12).index)
        subset = subset.sort_values("nes")
        colors = ["#2563EB" if value < 0 else "#D69E2E" for value in subset["nes"]]
        ax.barh(range(len(subset)), subset["nes"], color=colors)
        ax.set_yticks(range(len(subset)), subset["term"].map(shorten_term))
        ax.axvline(0, color="#374151", linewidth=0.8)
        ax.set_xlabel("Normalized enrichment score")
        ax.set_title(f"{contrast.replace('vsNC', '')}: top preranked GSEA results")

    ora = enrichment["ora"].copy()
    ora_sig = ora[ora["fdr"] < config["gsea"]["main_fdr"]].nsmallest(15, "fdr").copy()
    if ora_sig.empty:
        axes[1, 1].text(0.5, 0.5, "No ORA terms passed FDR < 0.05", ha="center", va="center")
        axes[1, 1].set_axis_off()
    else:
        ora_sig["score"] = -np.log10(ora_sig["fdr"].clip(lower=np.finfo(float).tiny))
        ora_sig["display_term"] = (
            ora_sig["contrast"].str.replace("vsNC", "", regex=False)
            + " "
            + ora_sig["direction"]
            + ": "
            + ora_sig["term"].map(shorten_term)
        )
        ora_sig = ora_sig.sort_values("score")
        y = np.arange(len(ora_sig))
        color_map = {"up": "#D69E2E", "down": "#2563EB"}
        axes[1, 1].scatter(
            ora_sig["score"],
            y,
            s=25 + 18 * ora_sig["overlap"],
            color=[color_map[value] for value in ora_sig["direction"]],
            alpha=0.85,
            edgecolor="white",
        )
        axes[1, 1].set_yticks(y, ora_sig["display_term"])
        axes[1, 1].set_xlabel("-log10 ORA FDR")
        axes[1, 1].set_title("Strongest directional DEG over-representation")
    fig.suptitle("Pathway-level treatment responses", fontsize=15, fontweight="bold")
    save_figure(fig, figure_dir / "04_pathway_enrichment")


def package_versions() -> dict[str, str]:
    packages = ["numpy", "pandas", "scipy", "scikit-learn", "statsmodels", "matplotlib", "seaborn", "gseapy"]
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def output_inventory(paths: Iterable[Path]) -> list[dict[str, Any]]:
    rows = []
    for path in sorted({candidate for candidate in paths if candidate.is_file()}):
        rows.append({"path": relative(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    return rows


def validate_report_links(report_dir: Path) -> list[Path]:
    reports = sorted(report_dir.glob("report_*.md"))
    image_pattern = re.compile(r"!\[[^\]]*\]\(([^)]+)\)")
    for report in reports:
        text = report.read_text(encoding="utf-8")
        for target in image_pattern.findall(text):
            clean_target = target.split("#", 1)[0]
            if clean_target.startswith(("http://", "https://")):
                continue
            candidate = (report.parent / clean_target).resolve()
            if not candidate.is_file():
                raise RuntimeError(f"Broken report image link in {report}: {target}")
    return reports


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True, help="Mentor-provided .tar.gz archive")
    parser.add_argument("--config", type=Path, default=ROOT / "analysis/config.json")
    parser.add_argument("--offline", action="store_true", help="Forbid reference downloads")
    parser.add_argument("--refresh-references", action="store_true", help="Replace cached reference snapshots")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    archive = args.archive.expanduser().resolve()
    config_path = args.config.expanduser().resolve()
    if not archive.is_file():
        raise FileNotFoundError(archive)
    config = load_config(config_path)

    processed_dir = ROOT / "data/processed"
    reference_dir = ROOT / "data/reference"
    table_dir = ROOT / "results/tables"
    figure_dir = ROOT / "results/figures"
    for directory in [processed_dir, reference_dir, table_dir, figure_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = [
        "FASTQ-level QC and alignment reports were not included in the supplied core-table archive.",
        "The mentor pipeline's final DEG annotation-combination step failed; core count and DE tables passed archive checks.",
        "Spi and Teg are coded labels because full compound identities were not supplied.",
        "No tumor-versus-normal HCC signature was supplied; disease reversal and therapeutic efficacy are out of scope.",
    ]
    started = now_iso()
    destination = ROOT / config["archive_destination"]
    print(f"[{now_iso()}] Verifying and preparing archive")
    dataset_dir, extraction = extract_archive(archive, destination)
    input_profile = validate_inputs(dataset_dir, config)
    idep_paths = prepare_idep_files(dataset_dir, processed_dir / "idep", config)
    write_preanalysis_inventory(dataset_dir, table_dir)

    print(f"[{now_iso()}] Preparing versioned gene mapping and pathway libraries")
    mapping, gmt_paths, reference_manifest = prepare_references(
        reference_dir, config, args.offline, args.refresh_references, warnings
    )

    print(f"[{now_iso()}] Running sample QC and differential-expression validation")
    qc = run_qc(dataset_dir, table_dir, config)
    deg_summary = reproduce_deg_tables(dataset_dir, config)
    deg_summary.to_csv(table_dir / "deg_summary.tsv", sep="\t", index=False)
    de_results = analyze_de(dataset_dir, mapping, processed_dir, table_dir, config)

    print(f"[{now_iso()}] Running local preranked GSEA and over-representation analysis")
    enrichment = run_enrichment(de_results, gmt_paths, table_dir, config)

    print(f"[{now_iso()}] Rendering report figures")
    plot_qc(qc, figure_dir)
    plot_deg(de_results, deg_summary, figure_dir, config)
    plot_concordance(de_results, figure_dir, config)
    plot_pathways(enrichment, figure_dir, config)
    report_paths = validate_report_links(ROOT / "reports")

    if qc["outliers"]["flagged_below_floor"].any():
        flagged = qc["outliers"].loc[qc["outliers"]["flagged_below_floor"], "sample"].tolist()
        warnings.append(f"Heuristic within-group correlation flag: {', '.join(flagged)}. No samples were removed.")

    output_paths = list(table_dir.glob("*.tsv")) + list(figure_dir.glob("*.png")) + list(figure_dir.glob("*.pdf")) + idep_paths
    output_paths += list((processed_dir / "ranked").glob("*.rnk"))
    output_paths += list((processed_dir / "differential_expression").glob("*.tsv"))
    output_paths += report_paths
    manifest = {
        "project": config["project"],
        "started_at": started,
        "completed_at": now_iso(),
        "status": "completed",
        "source_archive": {
            "path": str(archive),
            "bytes": archive.stat().st_size,
            "sha256": sha256_file(archive),
            **extraction,
        },
        "dataset_directory": relative(dataset_dir),
        "input_profile": input_profile,
        "parameters": config,
        "references": reference_manifest["references"],
        "mapping_minimum_coverage": float(de_results["mapping_summary"]["mapping_coverage"].min()),
        "software": {"python": sys.version.split()[0], "platform": platform.platform(), **package_versions()},
        "warnings": warnings,
        "outputs": output_inventory(output_paths),
    }
    manifest_path = ROOT / "results/run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[{now_iso()}] Completed successfully: {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
