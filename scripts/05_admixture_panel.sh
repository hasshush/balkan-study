#!/usr/bin/env bash
# Build the ADMIXTURE panel and run the K sweep (cached).
#
# Panel: the 70 GEO Balkan samples plus up to 20 reference samples per
# population from the QC'd merged dataset, LD-pruned at r^2 0.2.
# Runs: K = 2..7 with seed 42, plus K = 3 with seeds 1 and 2 to illustrate
# label switching.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/data/processed/nb03"
mkdir -p "$OUT"
ANALYSIS="$ROOT/data/processed/nb01/analysis"
cd "$ROOT"

if [[ ! -f "$OUT/panel.bed" ]]; then
  python - <<'PY'
import csv
from collections import defaultdict

groups = set(
    "Bosniac Bosnian_Serb Serbian Montenegrin Macedonian Kosovar Croatian Bulgarian "
    "Greek Albanian Romanian Moldavian Hungarian Gagauz Turkish TSI Sicilian Spanish "
    "French Basque Sardinian Ukrainian Belarusian Polish Russian Czech FIN CEU GBR IBS "
    "English Norwegian Icelandic Orcadian Lithuanian Estonian".split()
)
balkan = {"Bosniac", "Bosnian_Serb", "Serbian", "Montenegrin", "Macedonian", "Kosovar"}
rows = defaultdict(list)
for record in csv.DictReader(open("data/metadata/modern_merged_samples.tsv"), delimiter="\t"):
    if record["group"] in groups:
        rows[record["group"]].append(record["iid"])
keep = []
for group, iids in rows.items():
    limit = len(iids) if group in balkan else 20
    keep += [(group, iid) for iid in iids[:limit]]
with open("data/processed/nb03/keep.txt", "w") as handle:
    for group, iid in keep:
        handle.write(f"{group} {iid} {group}\n")
print(f"panel samples: {len(keep)}")
PY
  plink2 --bfile "$ANALYSIS" --keep "$OUT/keep.txt" --maf 0.05 \
    --indep-pairwise 50 5 0.3 --out "$OUT/prune" >/dev/null 2>&1
  plink2 --bfile "$ANALYSIS" --keep "$OUT/keep.txt" --extract "$OUT/prune.prune.in" \
    --make-bed --out "$OUT/panel" >/dev/null 2>&1
  echo "panel SNPs: $(wc -l < "$OUT/prune.prune.in")"
fi

run_admixture() {
  local k="$1" seed="$2"
  local prefix="$OUT/panel.${k}.s${seed}"
  if [[ -f "${prefix}.Q" ]]; then
    echo "cached: K=$k seed=$seed"
    return
  fi
  echo "running ADMIXTURE K=$k seed=$seed ..."
  (cd "$OUT" && admixture --cv -s "$seed" panel.bed "$k" > "log.${k}.s${seed}.txt" 2>&1)
  mv "$OUT/panel.${k}.Q" "${prefix}.Q"
  mv "$OUT/panel.${k}.P" "${prefix}.P"
}

for k in 2 3 4 5 6 7; do
  run_admixture "$k" 42
done
for seed in 1 2; do
  run_admixture 3 "$seed"
done
echo "done"
