#!/usr/bin/env python3
"""Extract ancient Balkan individuals from AADR 1240k into PLINK.

The 2023 Cell paper "A Genetic History of the Balkans from Roman Frontier to
Slavic Migrations" (Olalde, Carrion et al., doi:10.1016/j.cell.2023.10.018)
deposited raw sequencing data at ENA under PRJEB66422 and 1240k haploid
genotypes at the Reich lab. AADR v66.p1 incorporates those genotypes, so we
use the curated AADR copies here and cite the original study.

Selected individuals: any AADR group label from Croatia, Serbia, Montenegro,
North Macedonia, Albania, Slovenia or Bosnia-Herzegovina, excluding labels
that mark present-day samples ("Modern").

Genotypes use the tgeno format described in scripts/02_aadr_to_plink.py:
0 = homozygous ALT, 1 = het, 2 = homozygous REF, 3 = missing.
Output uses A1 = ALT, A2 = REF like the present-day panel.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import balkan_study_io as pio

ROOT = Path(__file__).resolve().parents[1]
AADR_DIR = ROOT / "data" / "reference" / "aadr"
GENO = AADR_DIR / "v66.p1_1240K.aadr.patch.PUB.geno"
SNP = AADR_DIR / "v66.p1_1240K.aadr.patch.PUB.snp"
IND = AADR_DIR / "v66.p1_1240K.aadr.patch.PUB.ind"
SITE_SOURCE = ROOT / "data" / "processed" / "nb01" / "analysis.bim"
OUT_PREFIX = ROOT / "data" / "ancient" / "aadr_ancient_balkans"
OUT_META = ROOT / "data" / "metadata" / "aadr_ancient_balkans_samples.tsv"

PREFIXES = ("Croatia_", "Serbia_", "Montenegro_", "NorthMacedonia_", "Albania_",
            "Slovenia_", "BosniaHerzegovina_", "Kosovo_")

TGENO_TO_PLINK = np.array([0, 2, 3, 1], dtype=np.uint8)


def main() -> None:
    OUT_PREFIX.parent.mkdir(parents=True, exist_ok=True)

    with open(GENO, "rb") as handle:
        header = handle.read(48).rstrip(b"\x00").decode("ascii").split()
        n_ind, n_snp = int(header[1]), int(header[2])
    print(f"1240k panel: {n_ind:,} individuals x {n_snp:,} sites")

    wanted = {line.split()[1] for line in open(SITE_SOURCE)}
    snp_records = []
    site_idx = []
    with open(SNP) as handle:
        for j, line in enumerate(handle):
            fields = line.split()
            if fields[0] in wanted:
                snp_records.append((fields[0], fields[1], int(fields[3]), fields[4].upper(), fields[5].upper()))
                site_idx.append(j)
    site_idx = np.array(site_idx)
    print(f"sites shared with the modern analysis dataset: {len(site_idx):,}")

    samples = []
    with open(IND) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 3:
                continue
            label = fields[2]
            if label.startswith(PREFIXES) and "Modern" not in label:
                samples.append({"iid": fields[0], "sex": fields[1], "label": label})
    print(f"ancient individuals selected: {len(samples):,}")

    bytes_per_ind = (n_snp + 3) // 4
    lut = np.array([[(b >> 6) & 3, (b >> 4) & 3, (b >> 2) & 3, b & 3] for b in range(256)],
                   dtype=np.uint8).reshape(256, 4)
    all_ids = [line.split()[0] for line in open(IND) if len(line.split()) >= 3]
    id_to_row = {name: i for i, name in enumerate(all_ids)}

    genotypes = np.empty((len(site_idx), len(samples)), dtype=np.uint8)
    with open(GENO, "rb") as handle:
        for col, sample in enumerate(samples):
            i = id_to_row[sample["iid"]]
            handle.seek(48 + i * bytes_per_ind)
            block = np.frombuffer(handle.read(bytes_per_ind), dtype=np.uint8)
            codes = lut[block].reshape(-1)[:n_snp]
            genotypes[:, col] = TGENO_TO_PLINK[codes[site_idx]]

    out_records = [
        {"chrom": chrom, "rsid": rsid, "cm": 0, "pos": pos, "a1": alt, "a2": ref}
        for rsid, chrom, pos, ref, alt in snp_records
    ]
    out_samples = [{"fid": s["label"], "iid": s["iid"], "sex": s["sex"]} for s in samples]

    pio.write_bed(str(OUT_PREFIX) + ".bed", genotypes)
    pio.write_bim(str(OUT_PREFIX) + ".bim", out_records)
    pio.write_fam(str(OUT_PREFIX) + ".fam", out_samples)

    with open(OUT_META, "w") as handle:
        handle.write("iid\tgroup\tsex\n")
        for sample in out_samples:
            handle.write(f"{sample['iid']}\t{sample['fid']}\t{sample['sex']}\n")

    by_label: dict[str, int] = {}
    for sample in out_samples:
        by_label[sample["fid"]] = by_label.get(sample["fid"], 0) + 1
    report = (f"sites: {len(site_idx):,}\nsamples: {len(samples):,}\n"
              f"missing rate: {(genotypes == 1).mean():.4f}\n\n"
              + "\n".join(f"{k}: {v}" for k, v in sorted(by_label.items())) + "\n")
    (ROOT / "data" / "processed" / "aadr_ancient_balkans_report.txt").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
