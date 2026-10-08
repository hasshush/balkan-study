#!/usr/bin/env bash
set -euo pipefail

echo "== poreklo environment =="
echo "root: $PWD"

for tool in plink2 admixture bcftools python jupyter; do
  if command -v "$tool" >/dev/null 2>&1; then
    printf '%-10s %s\n' "$tool" "$(command -v "$tool")"
  else
    printf '%-10s MISSING\n' "$tool"
  fi
done

echo
echo "== versions =="
plink2 --version 2>&1 | head -2 || true
admixture --version 2>&1 | grep -m1 'Version' || true
bcftools --version 2>&1 | head -2 || true
python - <<'PY' || true
import numpy, scipy, pandas, sklearn, matplotlib, seaborn, statsmodels
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("scikit-learn", sklearn.__version__)
print("matplotlib", matplotlib.__version__)
print("seaborn", seaborn.__version__)
print("statsmodels", statsmodels.__version__)
PY
