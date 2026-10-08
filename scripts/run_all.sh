#!/usr/bin/env bash
# End-to-end pipeline, everything cached. Run from inside `nix develop`.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

bash scripts/00_download_geo.sh
WITH_GENO=1 bash scripts/00_download_aadr.sh

python scripts/01_geo_to_plink.py
python scripts/02_aadr_to_plink.py
python scripts/03_merge_modern.py
python scripts/06_aadr_ancient_to_plink.py

python scripts/build_notebooks.py
bash scripts/05_admixture_panel.sh

for notebook in notebooks/*.ipynb; do
  echo "executing $notebook"
  jupyter nbconvert --to notebook --execute --inplace "$notebook"
done

python scripts/sanitize_notebooks.py
echo "all done"
