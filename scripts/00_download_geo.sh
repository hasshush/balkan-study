#!/usr/bin/env bash
# Download the raw GEO data used by scripts/01_geo_to_plink.py.
# Files land in data/raw/geo/ and are verified against known sha256 sums.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="$ROOT/data/raw/geo"
mkdir -p "$DEST"

fetch() {
  local url="$1" out="$2" sha="$3"
  if [[ -f "$DEST/$out" ]] && echo "$sha  $DEST/$out" | sha256sum -c - >/dev/null 2>&1; then
    echo "ok (cached): $out"
    return
  fi
  echo "downloading: $out"
  curl -sL --fail -o "$DEST/$out" "$url"
  echo "$sha  $DEST/$out" | sha256sum -c -
}

fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59032/suppl/GSE59032_processed.txt.gz" \
  "GSE59032_processed.txt.gz" \
  "dd03509282bc1ee5c0ed149f30e092195f4d8c822a9c8781149c70dd795038dc"

fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE59nnn/GSE59032/matrix/GSE59032_series_matrix.txt.gz" \
  "GSE59032_series_matrix.txt.gz" \
  "11d2184b4ddf85f33ed18c0de48aec682ac011859a123776c0235590d8eb49d3"

# Optional: full GPL8888 platform dump (2.7 GB, mostly unrelated samples).
# Only the platform table matters (ID, Intensity Only, SNP_ID) and it has no
# coordinates, so this file is not needed for the conversion.
if [[ "${WITH_GPL:-0}" == "1" ]]; then
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8888/soft/GPL8888_family.soft.gz" \
    "GPL8888_family.soft.gz" \
    "921edfab43fae9968d5a9233c636f3fe5b132326b4837275cebed0e3a6452f6e"
fi

echo "done"
