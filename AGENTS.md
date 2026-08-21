# DrugRep agent profile — draft

> Status: non-final and expected to evolve with the project.

## Project context

DrugRep investigates repurposing established drugs for hepatocellular carcinoma (HCC). Candidate ranking combines (1) a weighted average of six computational-pharmacology algorithms and (2) an AetherCell-derived anti-similarity score between virtual drug perturbations and the HCC disease state. Selected candidates are followed by experimental RNA-seq analysis: expression quality control, differential expression, and GSEA/GO/KEGG-style enrichment.

The current bioinformatics dataset contains four coded groups (`NC`, `Sora`, `Spi`, `Teg`) with three replicates per group. Do not expand ambiguous abbreviations into compound names without authoritative metadata. Do not claim disease reversal without a matched disease signature, or therapeutic efficacy from in-vitro transcriptomic association.

## Working rules

- Treat `data/raw/` as immutable. All transformations must be scripted and written to `data/processed/` or `results/`.
- Record input checksums, reference versions, parameters, seeds, software versions, and material warnings.
- Prefer full differential-expression tables for ranked analyses; never run GSEA only on thresholded DEGs.
- Use tested genes as the background for over-representation analysis.
- Separate verified results, interpretation, limitations, and unanswered mentor questions.
- Keep the workspace compact and update the relevant README when introducing a new durable directory or workflow.

## Important-task reports

For each important task, create a concise Markdown report unless the user explicitly waives it. The default filename is:

`reports/report_<YYYYMMDD-HHMM>_<snake_case_task_name>.md`

Use the Asia/Shanghai task timestamp. A user-specified name or format overrides this default. When needed, include the most decision-relevant plots, graphs, and small tables, each with a short caption or description. Keep the report focused; store exhaustive results in `results/tables/` and link them rather than pasting them into the report.

