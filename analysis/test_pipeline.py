import unittest

import pandas as pd

from analysis.run_all import build_ranked_table, normalize_gsea_result, run_ora, safe_tar_member


class PipelineUnitTests(unittest.TestCase):
    def test_archive_paths(self):
        self.assertTrue(safe_tar_member("dataset/table.tsv"))
        self.assertFalse(safe_tar_member("../escape.tsv"))
        self.assertFalse(safe_tar_member("/absolute/path.tsv"))

    def test_ranked_table_is_unique_and_directional(self):
        de = pd.DataFrame(
            {
                "stat": [4.0, -5.0, 1.0],
                "log2FoldChange": [1.2, -2.0, 0.2],
                "padj": [0.01, 0.001, 0.5],
            },
            index=["ENSG1", "ENSG2", "ENSG3"],
        )
        mapping = pd.DataFrame(
            {
                "ensembl_gene_id": ["ENSG1", "ENSG2", "ENSG3"],
                "gene_symbol": ["DUP", "DUP", "OTHER"],
            }
        )
        ranked = build_ranked_table(de, mapping)
        self.assertEqual(ranked["gene_symbol"].nunique(), 2)
        self.assertEqual(ranked.loc[ranked["gene_symbol"] == "DUP", "stat"].item(), -5.0)
        self.assertFalse(ranked["rank_score"].duplicated().any())

    def test_ora_uses_explicit_background(self):
        result = run_ora(
            {"A", "B", "C"},
            {"A", "B", "C", "D", "E", "F"},
            {"enriched": ["A", "B", "C"], "other": ["D", "E", "F"]},
            min_size=2,
            max_size=10,
        )
        self.assertEqual(result.iloc[0]["term"], "enriched")
        self.assertEqual(result.iloc[0]["background_size"], 6)
        self.assertGreaterEqual(result.iloc[0]["fdr"], result.iloc[0]["pvalue"])

    def test_gsea_result_is_canonicalized(self):
        raw = pd.DataFrame(
            {
                "Term": ["Z_PATHWAY", "A_PATHWAY"],
                "NES": [1.2, -1.5],
                "FDR q-val": [0.04, 0.01],
            },
            index=[17, 4],
        )
        normalized = normalize_gsea_result(raw)
        self.assertEqual(normalized["term"].tolist(), ["A_PATHWAY", "Z_PATHWAY"])
        self.assertNotIn("index", normalized.columns)


if __name__ == "__main__":
    unittest.main()
