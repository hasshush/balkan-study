# Provenance log

Record every download and tool version here.

## Tools (pinned in flake.lock, nixpkgs rev 7dd199b0e299)

| Tool | Version | Source |
|---|---|---|
| plink2 | 2.0.0-a.7.9 (27 Sep 2026) | nixpkgs |
| bcftools | 1.22 (htslib 1.22) | nixpkgs |
| admixture | 1.4.0 (static binary) | https://dalexander.github.io/admixture/binaries/admixture_linux-1.4.0.tar.gz, packaged in `nix/admixture.nix` |
| python | 3.14.7 (numpy 2.5.2, scipy 1.18.1, pandas 3.0.4, scikit-learn 1.9.0, statsmodels 0.14.6) | nixpkgs |
| jupyterlab | 4.6.4 | nixpkgs |

## GEO data (Kovacevic et al. 2014)

| Dataset | URL | sha256 | Notes |
|---|---|---|---|
| GSE59032 processed genotypes | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59032/suppl/GSE59032_processed.txt.gz | dd03509282bc1ee5c0ed149f30e092195f4d8c822a9c8781149c70dd795038dc | 809,245,435 bytes; GenomeStudio wide table; 70 samples |
| GSE59032 series matrix | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59032/matrix/GSE59032_series_matrix.txt.gz | 11d2184b4ddf85f33ed18c0de48aec682ac011859a123776c0235590d8eb49d3 | sample population/region metadata |
| GPL8888 platform soft | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8888/soft/GPL8888_family.soft.gz | 921edfab43fae9968d5a9233c636f3fe5b132326b4837275cebed0e3a6452f6e | probe IDs only, no coordinates; not required |

GEO sample composition (population/region fields): Bosniac 15, Bosnian Serb 4,
Serbian 14, Montenegrin 14, Macedonian 14, Kosovar 9. Croatian/Bosnian-Croat
samples from the paper are not in GEO.

## AADR v66.p1 (Reich lab, doi:10.7910/DVN/FFIDCW)

All downloaded files verified against the published `v66.p1__files.md5sum`.

| Dataset | Dataverse file id | Notes |
|---|---|---|
| v66.p1_1240K.aadr.patch.PUB.geno | 13994829 | 7,117,276,654 bytes, tgeno format |
| v66.p1_1240K.aadr.patch.PUB.snp | 13994514 | sha256 1189e253b82058822fddabc4966c1448dea70373a96c40106e720bedd93b55bf |
| v66.p1_1240K.aadr.patch.PUB.ind | 13994513 | |
| v66.p1_1240K.aadr.PUB.anno | 13994515 | |
| v66.p1_compatibility_HO.aadr.patch.PUB.snp | 13994524 | sha256 74582cfb2cde3943dcb969a3e80a2cdb711d6f2635e95c3adddaaf3a31fe1bbd |
| v66.p1_compatibility_HO.aadr.patch.PUB.ind | 13994523 | |
| v66.p1_compatibility_HO.aadr.PUB.anno | 13994525 | |
| v66.p1_HO.aadr.patch.PUB.geno | 13994808 | 4,029,634,650 bytes, tgeno format |
| v66.p1_HO.aadr.patch.PUB.snp | 13994527 | sha256 c0d565ee9e6a6edf30c2063430ee4f828e5db45f26d4747fd76f19eee650a5de |
| v66.p1_HO.aadr.patch.PUB.ind | 13994526 | |
| v66.p1_HO.aadr.PUB.anno | 13994528 | |
| v66.p1__files.md5sum | 13995422 | |
| aadr_v66.p1__README.docx | 13994530 | tgeno format description |

### tgeno format (validated)

ASCII 48-byte header `TGENO   <n_ind> <n_snp> <hash> <hash>` followed by
individual-major genotypes: for each individual, all SNPs packed 2 bits each,
4 SNPs per byte, first SNP of each group in the high bits. Codes:
`0 = homozygous ALT, 1 = het, 2 = homozygous REF, 3 = missing`.
Validated by 42/42 genotype concordance with 1000 Genomes TSI calls
(scripts/check_notes).

## Processed datasets

| Dataset | Script | Notes |
|---|---|---|
| data/processed/geo_balkan.{bed,bim,fam} | scripts/01_geo_to_plink.py | 530,659 SNPs, 70 GEO samples; alleles in arbitrary orientation |
| data/reference/ho_modern.{bed,bim,fam} | scripts/02_aadr_to_plink.py | 133,904 sites, 2,816 present-day samples, A1=ALT/A2=REF |
| data/processed/modern_merged.{bed,bim,fam} | scripts/03_merge_modern.py | 129,269 SNPs, 2,374 samples after removing 512 AADR technical duplicates; GEO oriented by frequency matching (r=0.967) |
| data/processed/nb01/analysis.{bed,bim,fam} | notebook 01 | QC'd analysis set: autosomes, MAF>1%, geno<5%, HWE 1e-6, LD-pruned; 72,552 SNPs x 2,374 samples |
| data/ancient/aadr_ancient_balkans.{bed,bim,fam} | scripts/06_aadr_ancient_to_plink.py | 679 ancient Balkan individuals at the 72,552 modern sites; A1=ALT/A2=REF |
| data/processed/nb03/panel.* | scripts/05_admixture_panel.sh | ADMIXTURE panel (509 samples, 59,278 SNPs); runs K=2..7 cached |

## External studies cited

| Study | Role | Accession / URL |
|---|---|---|
| Olalde, Carrion et al. 2023, Cell, "A Genetic History of the Balkans from Roman Frontier to Slavic Migrations" | ancient Balkan genomes (notebook 07) | ENA PRJEB66422; 1240k genotypes redistributed in AADR v66.p1; doi:10.1016/j.cell.2023.10.018 |
| Kovacevic et al. 2014, PLoS ONE | modern GEO samples | GEO GSE59032; doi:10.1371/journal.pone.0105090 |
| 1000 Genomes phase 3 | validation of tgeno decoding | EBI FTP release 20130502 |
| ADMIXTURE 1.4.0 | ancestry estimation | doi:10.1186/1471-2105-12-246 (1.3.0) |

