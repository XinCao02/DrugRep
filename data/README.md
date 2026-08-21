# Data

- `raw/`: verified, unpacked mentor inputs. These files are immutable and ignored by Git.
- `reference/`: cached Ensembl mapping and GMT gene-set snapshots. These are ignored by Git; provenance and checksums are copied into the run manifest.
- `processed/`: reproducible derivatives, including iDEP-ready files, mapped differential-expression tables, and ranked gene lists.

The source archive remains outside this directory and is never deleted or moved by the pipeline.

