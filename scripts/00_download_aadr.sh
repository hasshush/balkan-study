#!/usr/bin/env bash
# Download the AADR v66.p1 files used as reference panel and position map.
# Dataverse dataset: doi:10.7910/DVN/FFIDCW
#
# The .snp/.ind/.anno files are small and always fetched. The large .geno
# matrices are fetched only when WITH_GENO=1 (1240k by default).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$ROOT/data/reference/aadr"
mkdir -p "$DEST"

BASE="https://dataverse.harvard.edu/api/access/datafile"

fetch() {
  local id="$1" out="$2" sha="${3:-}"
  if [[ -f "$DEST/$out" ]]; then
    echo "cached: $out"
  else
    echo "downloading: $out"
    curl -sL --fail -o "$DEST/$out" "$BASE/$id"
  fi
  if [[ -n "$sha" ]]; then
    echo "$sha  $DEST/$out" | sha256sum -c -
  fi
}

fetch 13994514 "v66.p1_1240K.aadr.patch.PUB.snp" \
  "1189e253b82058822fddabc4966c1448dea70373a96c40106e720bedd93b55bf"
fetch 13994513 "v66.p1_1240K.aadr.patch.PUB.ind"
fetch 13994515 "v66.p1_1240K.aadr.PUB.anno"

fetch 13994524 "v66.p1_compatibility_HO.aadr.patch.PUB.snp" \
  "74582cfb2cde3943dcb969a3e80a2cdb711d6f2635e95c3adddaaf3a31fe1bbd"
fetch 13994523 "v66.p1_compatibility_HO.aadr.patch.PUB.ind"
fetch 13994525 "v66.p1_compatibility_HO.aadr.PUB.anno"

if [[ "${WITH_GENO:-0}" == "1" ]]; then
  fetch 13994829 "v66.p1_1240K.aadr.patch.PUB.geno"
fi

echo "done"
echo "sha256 of genotype matrices (record in provenance.md):"
[[ -f "$DEST/v66.p1_1240K.aadr.patch.PUB.geno" ]] && \
  sha256sum "$DEST/v66.p1_1240K.aadr.patch.PUB.geno"
