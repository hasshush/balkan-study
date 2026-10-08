#!/usr/bin/env python3
"""Merge the 70 GEO Balkan samples with the AADR present-day reference panel.

Why orientation is needed
-------------------------
The GEO processed genotypes report Illumina A/B allele classes, and we do not
have the 660W manifest, so when converting them we could only guess at which
DNA letter corresponds to class A versus class B (see
scripts/01_geo_to_plink.py). The AADR panel, by contrast, is stored with
A1 = ALT and A2 = REF.

For each SNP shared between the two datasets we compare the frequency of the
GEO A1 allele in the 70 Balkan samples with the frequency of the panel ALT
allele. If the GEO frequency matches the panel REF frequency instead, the GEO
alleles (and dosages) are flipped. This frequency-based orientation is the
standard remedy when array strand/allele information is unavailable.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import balkan_study_io as pio

ROOT = Path(__file__).resolve().parents[1]
GEO_PREFIX = ROOT / "data" / "processed" / "geo_balkan"
PANEL_PREFIX = ROOT / "data" / "reference" / "ho_modern"
OUT_PREFIX = ROOT / "data" / "processed" / "modern_merged"
OUT_META = ROOT / "data" / "metadata" / "modern_merged_samples.tsv"
REPORT = ROOT / "data" / "processed" / "merge_report.txt"

MIN_GEO_CALLS = 20
MIN_PANEL_CALLS = 50

# A1 dosage (-1 missing, 0, 1, 2) -> PLINK 2-bit code
DOSAGE_TO_CODE = pio.DOSAGE_TO_BED_CODE


def allele_set(record: dict) -> frozenset[str]:
    return frozenset({record["a1"], record["a2"]})


def main() -> None:
    geo_fam = pio.read_fam(str(GEO_PREFIX) + ".fam")
    geo_bim = pio.read_bim(str(GEO_PREFIX) + ".bim")
    geo_dos = pio.read_bed_dosage(str(GEO_PREFIX) + ".bed", len(geo_fam), len(geo_bim))
    print(f"GEO: {len(geo_bim):,} sites x {len(geo_fam)} samples")

    panel_fam = pio.read_fam(str(PANEL_PREFIX) + ".fam")
    panel_bim = pio.read_bim(str(PANEL_PREFIX) + ".bim")
    panel_dos = pio.read_bed_dosage(str(PANEL_PREFIX) + ".bed", len(panel_fam), len(panel_bim))
    print(f"panel: {len(panel_bim):,} sites x {len(panel_fam)} samples")

    # AADR lists the same individuals once per source assay (e.g. NA12878.DG
    # and NA12878.HO). Keep one copy per person, choosing the copy with the
    # fewest missing genotypes in this site set.
    def base_id(iid: str) -> str:
        return re.sub(r"\.(DG|SG|HO|AG|TW|BY|AA|EC|WGC)$", "", iid)

    panel_missing = (panel_dos < 0).sum(axis=0)
    groups: dict[str, list[int]] = {}
    for i, sample in enumerate(panel_fam):
        groups.setdefault(base_id(sample["iid"]), []).append(i)
    keep_mask = np.zeros(len(panel_fam), dtype=bool)
    n_duplicates = 0
    for indices in groups.values():
        best = min(indices, key=lambda i: int(panel_missing[i]))
        keep_mask[best] = True
        n_duplicates += len(indices) - 1
    panel_fam = [sample for sample, keep in zip(panel_fam, keep_mask) if keep]
    panel_dos = panel_dos[:, keep_mask]
    print(f"panel after de-duplication: {len(panel_fam)} samples ({n_duplicates} technical duplicates dropped)")

    geo_index = {record["rsid"]: i for i, record in enumerate(geo_bim)}

    kept_panel_rows: list[int] = []
    kept_geo_rows: list[int] = []
    n_no_match = n_allele_mismatch = n_low_calls = 0
    n_flips = 0

    for j, record in enumerate(panel_bim):
        k = geo_index.get(record["rsid"])
        if k is None:
            n_no_match += 1
            continue
        if allele_set(record) != allele_set(geo_bim[k]):
            n_allele_mismatch += 1
            continue
        geo_col = geo_dos[k]
        keep = geo_col >= 0
        if keep.sum() < MIN_GEO_CALLS:
            n_low_calls += 1
            continue
        panel_col = panel_dos[j]
        panel_keep = panel_col >= 0
        if panel_keep.sum() < MIN_PANEL_CALLS:
            n_low_calls += 1
            continue
        f_geo = geo_col[keep].mean() / 2.0
        f_panel_alt = panel_col[panel_keep].mean() / 2.0

        # The GEO A1 allele is either the panel ALT or the panel REF.
        if abs(f_geo - f_panel_alt) > abs((1.0 - f_geo) - f_panel_alt):
            geo_dos[k] = np.where(geo_col >= 0, 2 - geo_col, geo_col).astype(np.int8)
            n_flips += 1
        kept_panel_rows.append(j)
        kept_geo_rows.append(k)

    geo_rows = np.array(kept_geo_rows)
    panel_rows = np.array(kept_panel_rows)
    merged_dos = np.concatenate([geo_dos[geo_rows], panel_dos[panel_rows]], axis=1)
    merged_codes = DOSAGE_TO_CODE[(merged_dos + 1).astype(np.int8)]

    out_bim = [panel_bim[j] for j in kept_panel_rows]
    out_fam = (
        [{"fid": s["fid"], "iid": s["iid"], "sex": "0"} for s in geo_fam]
        + [{"fid": s["fid"], "iid": s["iid"], "sex": s.get("sex", "0")} for s in panel_fam]
    )

    pio.write_bed(str(OUT_PREFIX) + ".bed", merged_codes)
    pio.write_bim(str(OUT_PREFIX) + ".bim", out_bim)
    pio.write_fam(str(OUT_PREFIX) + ".fam", out_fam)

    OUT_META.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_META, "w") as handle:
        handle.write("iid\tgroup\tsource\tsex\n")
        for sample in geo_fam:
            handle.write(f"{sample['iid']}\t{sample['fid']}\tgeo\t0\n")
        for sample in panel_fam:
            handle.write(f"{sample['iid']}\t{sample['fid']}\taadr\t{sample.get('sex', '0')}\n")

    report = (
        f"panel sites:            {len(panel_bim):,}\n"
        f"not present in GEO:     {n_no_match:,}\n"
        f"allele set mismatches:  {n_allele_mismatch:,}\n"
        f"too few calls:          {n_low_calls:,}\n"
        f"sites kept:             {len(out_bim):,}\n"
        f"GEO orientation flips:  {n_flips:,} ({100.0 * n_flips / max(len(out_bim), 1):.1f}%)\n"
        f"duplicate individuals:  {n_duplicates:,}\n"
        f"samples in merge:       {len(out_fam):,}\n"
        f"missing rate:           {(merged_dos < 0).mean():.4f}\n"
    )
    REPORT.write_text(report)
    print(report)


if __name__ == "__main__":
    main()
