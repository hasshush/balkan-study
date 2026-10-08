#!/usr/bin/env python3
"""Convert GEO GSE59032 GenomeStudio genotypes to PLINK 1.x tped/tfam.

Background
----------
GSE59032 (Kovacevic et al. 2014) contains 70 western Balkan individuals
genotyped on the Illumina Human660W-Quad v1.0 array. GEO stores the output
of Illumina GenomeStudio as a wide table: one row per probe, and six columns
per sample (GType, Score, Theta, R, B Allele Freq, Log R Ratio).

GType is a two-letter genotype call built from the array's two alleles:
    AA = two copies of allele A, AB = one A and one B, BB = two B.

The probe identifiers (`ID_REF`) are rsIDs for SNP probes (and `cnv...` for
copy-number probes, which are skipped). GEO's platform annotation (GPL8888)
does not provide coordinates, so genomic positions and alleles are taken
from the AADR v66.p1 .snp maps (GRCh37). Probes without an rsID in any AADR
map are skipped.

Important: the array's A/B labels are not guaranteed to correspond to the
AADR REF/ALT alleles in any particular order. This script writes REF first
and ALT second; orientation is corrected later when merging with the
reference panel (see scripts/03_merge_reference.py and the notebooks).
Within this dataset alone, allele orientation is irrelevant for
ADMIXTURE/PCA/FST-style analyses.

Outputs
-------
data/processed/geo_balkan.tped   PLINK transposed text genotypes
data/processed/geo_balkan.tfam   family/population labels
data/processed/geo_conversion_report.txt
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEO_PROCESSED = ROOT / "data" / "raw" / "geo" / "GSE59032_processed.txt.gz"
GEO_MATRIX = ROOT / "data" / "raw" / "geo" / "GSE59032_series_matrix.txt.gz"

# Order matters: first map wins for a given rsID.
AADR_SNP_MAPS = [
    ROOT / "data" / "reference" / "aadr" / "v66.p1_1240K.aadr.patch.PUB.snp",
    ROOT / "data" / "reference" / "aadr" / "v66.p1_compatibility_HO.aadr.patch.PUB.snp",
    ROOT / "data" / "reference" / "aadr" / "v66.p1_HO.aadr.patch.PUB.snp",
]

OUT_TPED = ROOT / "data" / "processed" / "geo_balkan.tped"
OUT_TFAM = ROOT / "data" / "processed" / "geo_balkan.tfam"
OUT_REPORT = ROOT / "data" / "processed" / "geo_conversion_report.txt"


def load_snp_map() -> dict[str, tuple[str, str, str, str]]:
    """rsID -> (chrom, pos, ref, alt), GRCh37, from AADR .snp files."""
    snp_map: dict[str, tuple[str, str, str, str]] = {}
    for path in AADR_SNP_MAPS:
        if not path.exists():
            print(f"warning: missing {path}", file=sys.stderr)
            continue
        added = 0
        with open(path) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 6:
                    continue
                rsid, chrom, _genpos, pos, ref, alt = fields[:6]
                if rsid not in snp_map:
                    snp_map[rsid] = (chrom, pos, ref, alt)
                    added += 1
        print(f"loaded {added:>9,} new markers from {path.name}")
    return snp_map


def load_sample_metadata() -> list[tuple[str, str, str]]:
    """Parse the series matrix into (title, population, region) in column order."""
    titles: list[str] = []
    populations: list[str] = []
    regions: list[str] = []
    with gzip.open(GEO_MATRIX, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                titles = [value.strip('"') for value in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                values = [value.strip('"') for value in line.rstrip("\n").split("\t")[1:]]
                if values and values[0].startswith("population:"):
                    populations = [value.split(":", 1)[1].strip() for value in values]
                elif values and values[0].startswith("region:"):
                    regions = [value.split(":", 1)[1].strip() for value in values]
    if not (titles and populations and regions):
        raise RuntimeError("could not parse sample metadata from series matrix")
    return list(zip(titles, populations, regions))


def group_label(population: str, region: str) -> str:
    """Collapse GEO population/region fields into analysis groups."""
    if population == "Serbian":
        if "Bosnia" in region:
            return "Bosnian_Serb"
        return "Serbian"
    return population


GENOTYPE_TOKENS = {
    "AA": ("ref", "ref"),
    "AB": ("ref", "alt"),
    "BB": ("alt", "alt"),
}


def main() -> None:
    snp_map = load_snp_map()
    metadata = load_sample_metadata()
    print(f"{len(metadata)} samples in series matrix")

    OUT_TPED.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(GEO_PROCESSED, "rt", errors="replace") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        gtype_columns = [i for i, name in enumerate(header) if name.endswith(".GType")]
        sample_names = [header[i][: -len(".GType")] for i in gtype_columns]
        if len(sample_names) != len(metadata):
            raise RuntimeError(
                f"sample count mismatch: {len(sample_names)} genotype columns "
                f"vs {len(metadata)} metadata rows"
            )
        if sample_names != [title for title, _, _ in metadata]:
            print("warning: genotype column order differs from series matrix order", file=sys.stderr)

        with open(OUT_TFAM, "w") as tfam, open(OUT_TPED, "w") as tped:
            for title, population, region in metadata:
                tfam.write(f"{group_label(population, region)}\t{title}\t0\t0\t0\t-9\n")

            n_rows = n_written = n_unmapped = n_missing_cells = n_duplicates = 0
            seen: set[str] = set()

            for line in handle:
                fields = line.split("\t")
                rsid = fields[0]
                n_rows += 1
                if not rsid.startswith("rs") or rsid not in snp_map:
                    n_unmapped += 1
                    continue
                if rsid in seen:
                    n_duplicates += 1
                    continue
                seen.add(rsid)

                chrom, pos, ref, alt = snp_map[rsid]
                tokens: list[str] = []
                ref_ref = (ref, ref)
                ref_alt = (ref, alt) if ref < alt else (alt, ref)
                alt_alt = (alt, alt)
                for i in gtype_columns:
                    genotype = fields[i]
                    if genotype == "AA":
                        pair = ref_ref
                    elif genotype == "AB":
                        pair = ref_alt
                    elif genotype == "BB":
                        pair = alt_alt
                    else:
                        pair = ("0", "0")
                        n_missing_cells += 1
                    tokens.extend(pair)
                tped.write(f"{chrom}\t{rsid}\t0\t{pos}\t" + "\t".join(tokens) + "\n")
                n_written += 1

    total_cells = n_written * len(sample_names)
    report = (
        f"probe rows read:            {n_rows:,}\n"
        f"markers written:            {n_written:,}\n"
        f"rows skipped (no AADR map): {n_unmapped:,}\n"
        f"duplicate rsIDs skipped:    {n_duplicates:,}\n"
        f"samples:                    {len(sample_names)}\n"
        f"genotype cells:             {total_cells:,}\n"
        f"missing cells:              {n_missing_cells:,} "
        f"({100.0 * n_missing_cells / max(total_cells, 1):.2f}%)\n"
    )
    OUT_REPORT.write_text(report)
    print(report)
    print(f"wrote {OUT_TPED} and {OUT_TFAM}")


if __name__ == "__main__":
    main()
