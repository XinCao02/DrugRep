# Analysis

`run_all.py` is the single entry point for safe extraction, validation, quality control, differential-expression robustness checks, local enrichment, plotting, and manifest generation. Parameters are stored in `config.json`.

```bash
python analysis/run_all.py --archive /path/to/archive.tar.gz --config analysis/config.json
```

Use `--offline` only after references have been cached. Use `--refresh-references` to intentionally replace the cached mapping and GMT snapshots. The script is idempotent: verified inputs and reference files are reused.

