# poreklo

Genetic structure of Western Balkan populations: do self-identified nations
correspond to genetic groupings?

The project is a learning-oriented walkthrough that builds up to a test of
whether nation membership can be predicted from autosomal genotype data
better than chance, once geography is accounted for.

## Research question

- How is genetic variation structured among Western Balkan populations?
- Do self-identified nations (Bosniacs, Bosnian Serbs, Serbians,
  Montenegrins, Macedonians, Kosovars, Croats, ...) correspond to genetic
  clusters?
- Can nation be predicted from genotype, and does nation explain variation
  beyond geographic distance?

## Data

| Dataset | Role | Source |
|---|---|---|
| Kovacevic et al. 2014, GSE59032 | 70 western Balkan individuals, Illumina Human660W-Quad | NCBI GEO |
| AADR (Reich lab) | modern reference populations (Europe/Balkans) | Harvard Dataverse |
| Carrión et al. 2023 ancients (PMC10752003) | ancient Balkan genomes, kept separate | ENA |

Details, hashes and versions are recorded in `provenance.md`.

Note: GEO GSE59032 contains only the 70 newly genotyped samples. The 20
Croatians/Bosnian Croats and the global reference samples used in the paper
come from other published sources.

## Environment

NixOS with flakes enabled. From the repository root:

```
nix develop          # or: nix-shell
bash scripts/check_env.sh
```

The shell provides `plink2`, `admixture` (packaged in `nix/admixture.nix`),
`bcftools`, and a Python environment with numpy/scipy/pandas/scikit-learn/
matplotlib/seaborn/statsmodels/jupyterlab.

## Layout

```
nix/            package definitions (admixture)
scripts/        download, conversion and QC scripts
notebooks/      numbered walkthrough (01 .. 07)
data/raw/       untouched downloads (gitignored)
data/processed/ PLINK datasets and metadata (gitignored)
data/reference/ AADR / reference panels (gitignored)
data/ancient/   2023 Cell ancient data, analyzed separately (gitignored)
results/        figures and tables (gitignored)
provenance.md   data provenance log
```

## Walkthrough

1. `01_data_foundations` - genotype arrays, dosages, QC, conversion to PLINK
2. `02_pca` - PCA intuition via simulations, then real data
3. `03_admixture` - model-based soft clustering, K selection
4. `04_fst_amova` - differentiation between populations
5. `05_kmeans` - k-means on PCs, K selection, permutation nulls
6. `06_nation_test` - supervised prediction and geography controls
7. `07_ancient` - separate track for ancient samples

## Working notes

- Ancient samples are never merged into the modern analysis.
- All random steps use recorded seeds.
- Balkan FST values are tiny (~0.0005-0.001); results are reported with
  permutation tests and confidence intervals, not just point estimates.
