# Provenance log

Record every download and tool version here.

## Tools

| Tool | Version | Source |
|---|---|---|
| plink2 | reported by `plink2 --version` | nixpkgs |
| bcftools | reported by `bcftools --version` | nixpkgs |
| admixture | 1.4.0 (static binary) | https://dalexander.github.io/admixture/binaries/admixture_linux-1.4.0.tar.gz, packaged in `nix/admixture.nix` |
| Python | nixpkgs `python3` | nixpkgs |

## Downloads

| Date | Dataset | URL | sha256 | Notes |
|---|---|---|---|---|
| 2026-10-07 | GSE59032 processed genotypes | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59032/suppl/GSE59032_processed.txt.gz | TODO | 809,245,435 bytes, GenomeStudio wide table |
| 2026-10-07 | GSE59032 signal intensities | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59032/suppl/GSE59032_signal_intensities.txt.gz | TODO | 598,584,810 bytes |
| 2026-10-07 | GSE59032 series matrix | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59032/matrix/GSE59032_series_matrix.txt.gz | TODO | sample metadata |
| 2026-10-07 | GPL8888 platform | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8888/soft/GPL8888_family.soft.gz | TODO | probe IDs only; no coordinates |
| 2026-10-07 | ADMIXTURE 1.4.0 linux binary | https://dalexander.github.io/admixture/binaries/admixture_linux-1.4.0.tar.gz | sha256-VSDE8YUl9JqyxxXXEbQPSggY7xSFNhqU7IlU0TZObxk= | statically linked |

## Pending decisions

- Illumina Human660W-Quad v1.0 manifest source (for probe -> rsID/chr/pos).
- AADR release/version to use.
- 2023 Cell ancient data release files to download.
