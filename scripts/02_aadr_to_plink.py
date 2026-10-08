#!/usr/bin/env python3
"""Convert AADR 'tgeno' genotypes to PLINK for a curated present-day panel.

AADR v66.p1 stores genotypes in the 'transpose_packed' (tgeno) format:

    bytes 0..47   ASCII header: "TGENO   <n_ind> <n_snp> <hash> <hash>" + padding
    after that    individual-major genotypes: for each individual, the genotypes
                  for all SNPs packed two bits each, four SNPs per byte, with the
                  first SNP of each group in the two HIGH bits of the byte

Codes are 0 = homozygous ALT, 1 = heterozygous, 2 = homozygous REF, 3 = missing
(validated against 1000 Genomes TSI genotypes; note this is the reverse of the
classic EIGENSTRAT convention).

This script selects present-day populations and writes a PLINK 1 dataset with
A1 = ALT (variant) and A2 = REF, restricted to SNPs shared with the GEO data.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import poreklo_io as pio

ROOT = Path(__file__).resolve().parents[1]
AADR_DIR = ROOT / "data" / "reference" / "aadr"
GENO = AADR_DIR / "v66.p1_HO.aadr.patch.PUB.geno"
SNP = AADR_DIR / "v66.p1_HO.aadr.patch.PUB.snp"
IND = AADR_DIR / "v66.p1_HO.aadr.patch.PUB.ind"
GEO_RSIDS = ROOT / "data" / "processed" / "geo_rsids.txt"
OUT_PREFIX = ROOT / "data" / "reference" / "ho_modern"
OUT_META = ROOT / "data" / "metadata" / "ho_modern_samples.tsv"

# Present-day populations to keep for the reference panel. Labels are exact
# AADR .ind group labels; ancient samples use country_period style labels and
# are excluded by construction.
POPULATIONS = [
    # Western Balkans and immediate neighbours
    "Albanian", "Bulgarian", "Croatian", "Greek", "Romanian", "Moldavian",
    "Hungarian", "Gagauz", "Turkish", "Cypriot",
    # Southern Europe
    "Spanish", "French", "Basque", "Sardinian", "Sicilian", "TSI", "IBS", "Maltese",
    # Central, Eastern and Northern Europe
    "English", "Scottish", "Orcadian", "Icelandic", "Norwegian", "Czech", "Polish",
    "Ukrainian", "Belarusian", "Russian", "Lithuanian", "Estonian", "Finnish",
    "FIN", "CEU", "GBR",
    # Caucasus
    "Abkhasian", "Georgian", "Armenian", "Chechen", "Adygei", "Lezgin", "Ossetian",
    "Balkar", "Kumyk", "Darginian", "Kabardinian", "Ingushian", "Avar",
    "Circassian", "Abazin", "Karachai", "Tabasaran", "Lak", "Kaitag", "Kubachinian",
    # Middle East
    "Iranian", "Assyrian", "Druze", "Palestinian", "Jordanian", "Lebanese",
    "Syrian", "Egyptian", "Saudi", "Yemeni", "Emirati", "BedouinA", "BedouinB",
    "Samaritan",
    # global outgroups, for context in the global PCA
    "YRI", "Han", "Japanese", "CHB", "JPT", "CHS", "Papuan", "Mbuti", "Biaka",
    "Karitiana", "Pima", "Mayan", "Surui",
]

# tgeno code -> PLINK 2-bit code, given A1 = ALT, A2 = REF
# tgeno 0 (hom ALT) -> 00 (hom A1); 1 (het) -> 10; 2 (hom REF) -> 11; 3 (missing) -> 01
TGENO_TO_PLINK = np.array([0, 2, 3, 1], dtype=np.uint8)


def parse_header(handle) -> int:
    raw = handle.read(48)
    if raw[:5] != b"TGENO":
        raise RuntimeError("not a tgeno file")
    text = raw.rstrip(b"\x00").decode("ascii")
    fields = text.split()
    n_snp = int(fields[2])
    print(f"tgeno header: n_ind={fields[1]} n_snp={n_snp} hashes={fields[3:5]}")
    return n_snp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-sites", action="store_true", help="keep all HO sites instead of only GEO-shared ones")
    args = parser.parse_args()

    n_snp = None
    with open(GENO, "rb") as handle:
        n_snp = parse_header(handle)

    snp_records = []
    with open(SNP) as handle:
        for line in handle:
            fields = line.split()
            snp_records.append((fields[0], fields[1], int(fields[3]), fields[4].upper(), fields[5].upper()))
    if len(snp_records) != n_snp:
        raise RuntimeError(f"snp count mismatch: {len(snp_records)} vs {n_snp}")

    samples = []
    with open(IND) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 3:
                samples.append({"iid": fields[0], "sex": fields[1], "label": fields[2]})
    if len(samples) != int(open(GENO, "rb").read(48).rstrip(b"\x00").decode().split()[1]):
        raise RuntimeError("ind count does not match tgeno header")

    preferred = set(open(GEO_RSIDS).read().split()) if not args.all_sites else None
    site_mask = np.zeros(n_snp, dtype=bool)
    for j, (rsid, _chrom, _pos, ref, alt) in enumerate(snp_records):
        if len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT":
            if preferred is None or rsid in preferred:
                site_mask[j] = True
    site_idx = np.flatnonzero(site_mask)
    print(f"sites selected: {len(site_idx):,}")

    selected = [i for i, s in enumerate(samples) if s["label"] in set(POPULATIONS)]
    print(f"individuals selected: {len(selected):,}")

    bytes_per_ind = (n_snp + 3) // 4
    lut = np.array([[(b >> 6) & 3, (b >> 4) & 3, (b >> 2) & 3, b & 3] for b in range(256)], dtype=np.uint8).reshape(256, 4)
    genotypes = np.empty((len(site_idx), len(selected)), dtype=np.uint8)
    with open(GENO, "rb") as handle:
        for col, i in enumerate(selected):
            handle.seek(48 + i * bytes_per_ind)
            block = np.frombuffer(handle.read(bytes_per_ind), dtype=np.uint8)
            if block.size != bytes_per_ind:
                raise RuntimeError(f"short read for individual {i}")
            codes = lut[block].reshape(-1)[:n_snp]
            genotypes[:, col] = TGENO_TO_PLINK[codes[site_idx]]

    out_records = [
        {
            "chrom": snp_records[j][1],
            "rsid": snp_records[j][0],
            "cm": 0,
            "pos": snp_records[j][2],
            "a1": snp_records[j][4],  # ALT
            "a2": snp_records[j][3],  # REF
        }
        for j in site_idx
    ]
    out_samples = [
        {"fid": samples[i]["label"], "iid": samples[i]["iid"], "sex": samples[i]["sex"]}
        for i in selected
    ]

    pio.write_bed(str(OUT_PREFIX) + ".bed", genotypes.astype(np.uint8))
    pio.write_bim(str(OUT_PREFIX) + ".bim", out_records)
    pio.write_fam(str(OUT_PREFIX) + ".fam", out_samples)

    OUT_META.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_META, "w") as handle:
        handle.write("iid\tlabel\tsex\n")
        for sample in out_samples:
            handle.write(f"{sample['iid']}\t{sample['fid']}\t{sample['sex']}\n")

    missing = (genotypes == 1).mean()
    by_label: dict[str, int] = {}
    for sample in out_samples:
        by_label[sample["fid"]] = by_label.get(sample["fid"], 0) + 1
    report = (
        f"sites: {len(site_idx):,}\n"
        f"samples: {len(selected):,}\n"
        f"missing rate: {missing:.4f}\n\n"
        + "\n".join(f"{k}: {v}" for k, v in sorted(by_label.items()))
        + "\n"
    )
    (ROOT / "data" / "processed" / "ho_modern_report.txt").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
