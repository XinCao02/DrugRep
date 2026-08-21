# GitHub initial publish

## Outcome

Prepared the complete DrugRep workspace for publication to `xincao02/DrugRep` on the feature branch `codex/rnaseq-workspace`, preserving the existing remote `main` history.

## Included

- Reproducible RNA-seq analysis code, configuration, compact processed inputs, result tables, figures, and manifests.
- Project documentation, technical reports, the AetherCell reference paper, PowerPoint deck, animated HTML, and standalone offline HTML.
- Regeneration scripts for the presentation and standalone HTML.

## Data policy

The mentor archive and unpacked raw data, cached reference files, `.venv`, Python caches, and transient logs remain excluded by `.gitignore`. No tracked file exceeds 3 MB.

## Validation

- `python -m unittest analysis.test_pipeline -v`: 4/4 tests passed.
- Standalone presentation regenerated with all five images embedded.
- Remote repository inspected before integration; its existing `main` branch contained only `README.md`.

## Notes

The commit author is configured locally as `xincao02 <xincao02@users.noreply.github.com>`. This report records preparation; the final commit hash and push status are reported in the task handoff.
