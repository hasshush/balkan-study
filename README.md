# poreklo

Do self-identified nations in the Western Balkans correspond to genetic
groupings? This project is a learning-oriented walkthrough that builds from raw
genotype data to a direct statistical test of that question.

## Data

| Dataset | Role | Source |
|---|---|---|
| Kovacevic et al. 2014 (GSE59032) | 70 western Balkan individuals: Bosniac, Bosnian Serb, Serbian, Montenegrin, Macedonian, Kosovar | NCBI GEO |
| AADR v66.p1 (Reich lab) | 2,304 present-day reference samples: Croatian, Bulgarian, Greek, Albanian, Romanian, Moldavian, Hungarian, Gagauz, Turkish and wider Europe / West Asia | Harvard Dataverse, doi:10.7910/DVN/FFIDCW |
| AADR v66.p1 + Olalde & Carrion et al. 2023 | 679 ancient Balkan individuals (1st millennium CE), kept in a **separate** track | ENA PRJEB66422 / AADR |

The merged modern dataset contains **129,269 SNPs x 2,374 samples** (70 GEO +
2,304 reference). Details and checksums: `provenance.md`.

## Key findings

- **No separable nation clusters.** PCA shows one continuous north-west to
  south-east cline shared by all six nations; Kosovars sit at one edge,
  Bosniacs/Croats at the other. ADMIXTURE components are shared across nations.
- **Nation labels carry a weak but real signal.** Between-nation Hudson FST is
  ~0.0014 (permutation p ~ 0.001); random splits *within* a nation give FST
  values near zero. This is 1-2 orders of magnitude below between-country
  differentiation elsewhere in Europe.
- **The signal is weak, not decisive.** K-means partitions align with nation
  labels only slightly (ARI 0.04-0.13, p <= 0.03 depending on K).
- **Classification is above chance but far from perfect.** Logistic regression /
  random forest on PCA scores reach balanced accuracy ~0.45-0.47 against a
  chance level of 0.167 (p ~ 0.001). The same pipeline separates six reference
  populations at ~0.96. People are closer on average to their own nation's
  samples (p ~ 0.02), but the effect is tiny.
- **Mostly geography.** Sampling towns differ between nations, so the weak
  predictability is entangled with place; the Mantel tests with coarse
  group-level coordinates cannot separate the two.
- **One shared gene pool.** Genetics neither validates nor invalidates national
  identities: all six groups sit on the same South Slavic cline, which
  ancient DNA shows was shaped by local continuity plus Roman-era and
  Migration-period gene flow.

## Notebooks

Run `nix develop` (or `nix-shell`) from the repo root, then `jupyter lab`.

| Notebook | Contents |
|---|---|
| `01_data_foundations` | genotype arrays, dosages, missingness, MAF, HWE, LD pruning, relatedness, heterozygosity |
| `02_pca` | PCA intuition with simulations (clusters vs cline), real PCA of the Balkans |
| `03_admixture` | soft clustering, K selection by CV, label switching |
| `04_fst_amova` | Hudson FST, permutation test of nation labels, within vs between |
| `05_kmeans` | your original idea: k-means on PCs, K selection, ARI vs labels, permutation nulls |
| `06_nation_test` | supervised classification with cross-validation, own-group distance test, Mantel / partial Mantel |
| `07_ancient` | **separate track**: 679 ancient Balkan genomes projected next to the modern cline |

## Reproduce

```bash
nix develop
bash scripts/00_download_geo.sh                 # GEO genotypes + metadata (~820 MB)
WITH_GENO=1 bash scripts/00_download_aadr.sh    # AADR .snp/.ind/.anno/.geno (~11 GB)
python scripts/01_geo_to_plink.py               # GenomeStudio -> PLINK (530,659 SNPs)
python scripts/02_aadr_to_plink.py              # AADR tgeno -> PLINK panel
python scripts/03_merge_modern.py               # frequency-based orientation + merge
python scripts/06_aadr_ancient_to_plink.py      # ancient track
python scripts/build_notebooks.py               # notebooks_src/*.py -> notebooks/*.ipynb
bash scripts/05_admixture_panel.sh              # cached ADMIXTURE K=2..7
jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb
```

Notes:

- `scripts/01_geo_to_plink.py` handles GEO's GenomeStudio-format output; the
  AADR `tgeno` format is decoded in `scripts/02_aadr_to_plink.py` (format
  validated against 1000 Genomes genotypes, see `provenance.md`).
- GEO allele orientation is unknown without the Illumina manifest; it is
  resolved by matching allele frequencies against the reference panel
  (`scripts/03_merge_modern.py`, frequency correlation r = 0.967).
- All scripts cache their outputs; re-running is cheap after the first pass.

## Layout

```
nix/admixture.nix       packaged ADMIXTURE 1.4.0
scripts/                download, format conversion, merge, notebook builder
notebooks_src/          notebook sources (plain Python, # %% cells)
notebooks/              generated and executed notebooks
data/raw|reference|processed|ancient|metadata/
results/figures|tables/
provenance.md           downloads, hashes, format notes
```
