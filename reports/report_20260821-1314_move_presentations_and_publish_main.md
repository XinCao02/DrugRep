# Move presentations under reports and publish main

## Outcome

Moved the complete advisor-facing presentation workspace from `presentations/` to `reports/presentations/`. Updated repository documentation, report links, local-server URLs, and regeneration commands to the new location.

## Reproducibility fix

After the move, `build_rnaseq_deck.py` now resolves the project root two levels above the script while writing outputs beside the script. This preserves access to `results/tables/` and keeps generated presentation assets under `reports/presentations/`.

## Validation

- Rebuilt the six-slide PowerPoint and all presentation assets successfully.
- Regenerated the standalone HTML with five embedded images.
- RNA-seq pipeline unit tests: 4/4 passed.
- Confirmed that remote `main` already contained the prior project merge before preparing this change.

The final commit and remote verification are recorded in the task handoff.
