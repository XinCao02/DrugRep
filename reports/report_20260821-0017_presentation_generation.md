# RNA-seq advisor presentation

Created a six-slide, primarily Chinese PowerPoint for discussion with ZZH:

[`HCC_RNAseq_downstream_advisor_20260821.pptx`](presentations/HCC_RNAseq_downstream_advisor_20260821.pptx)

## Slide structure

1. Project title and interpretation boundary.
2. Outline plus the completed downstream workflow: input checks, expression QC, DEG analysis, robustness checks, GSEA/ORA, and cross-treatment comparison.
3. Raw-data inventory and QC, including what ZZH completed, what was independently verified, and which read-level artifacts remain unavailable.
4. RSEM–DESeq2 DEG totals and three volcano plots.
5. Cross-method robustness and candidate similarity to sorafenib.
6. Selected Hallmark GSEA results, ORA confirmation, conclusions, limitations, and questions for ZZH.

The deck uses a restrained academic style with an ivory background, navy/teal typography, and gold accents. All five charts were regenerated from validated project TSV files rather than copied from the Markdown report.

## Validation

- Six slides, satisfying the requested limit of fewer than seven.
- Chinese is the primary language; all references to the mentor use `ZZH`.
- The PPTX ZIP container, embedded media, slide bounds, and text-density guardrail are checked automatically.
- The deck explicitly states that treatment similarity is descriptive and that no tumor-versus-normal signature was supplied, so HCC reversal and efficacy are not established.

Rebuild with:

```bash
.venv/bin/python reports/presentations/build_rnaseq_deck.py
```
