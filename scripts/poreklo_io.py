#!/usr/bin/env python3
"""Small numpy helpers for PLINK 1 binary files.

Facts about the .bed format (PLINK 1):
- the first three bytes are the magic number 0x6c 0x1b 0x01
- after that, genotypes are stored SNP-major: for each SNP, all samples in order,
  packed four samples per byte
- within a byte, the first sample occupies the two lowest bits
- the 2-bit codes are:
      00 -> homozygous for allele A1 (dosage 2 of A1)
      01 -> missing
      10 -> heterozygous (dosage 1 of A1)
      11 -> homozygous for allele A2 (dosage 0 of A1)

The helper functions here map between those codes and additive dosages of the
A1 allele, with -1 for missing.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

PLINK_MAGIC = b"\x6c\x1b\x01"

# 2-bit .bed code -> A1 dosage; code 1 (missing) -> -1
BED_CODE_TO_DOSAGE = np.array([2, -1, 1, 0], dtype=np.int8)

# A1 dosage -> 2-bit .bed code, indexed by dosage + 1
# dosage -1 (missing) -> 1, 0 -> 3, 1 -> 2, 2 -> 0
DOSAGE_TO_BED_CODE = np.array([1, 3, 2, 0], dtype=np.uint8)


def read_fam(path: str | Path) -> list[dict]:
    """Read a .fam file into dicts with fid, iid and sex."""
    samples = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 2:
                continue
            samples.append(
                {
                    "fid": fields[0],
                    "iid": fields[1],
                    "sex": fields[4] if len(fields) > 4 else "0",
                }
            )
    return samples


def read_bim(path: str | Path) -> list[dict]:
    """Read a .bim file into dicts with chrom, rsid, pos, a1, a2."""
    records = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 6:
                continue
            records.append(
                {
                    "chrom": fields[0],
                    "rsid": fields[1],
                    "cm": fields[2],
                    "pos": int(fields[3]),
                    "a1": fields[4],
                    "a2": fields[5],
                }
            )
    return records


def read_bed_dosage(path: str | Path, n_samples: int, n_sites: int) -> np.ndarray:
    """Read a .bed file and return A1 dosages as an int8 matrix.

    Shape is (n_sites, n_samples); values are 0, 1, 2, with -1 for missing.
    """
    bytes_per_snp = (n_samples + 3) // 4
    raw = np.fromfile(path, dtype=np.uint8, offset=len(PLINK_MAGIC))
    expected = n_sites * bytes_per_snp
    if raw.size != expected:
        raise ValueError(f"{path}: expected {expected} bytes of genotypes, found {raw.size}")
    packed = raw.reshape(n_sites, bytes_per_snp)
    codes = np.empty((n_sites, bytes_per_snp * 4), dtype=np.uint8)
    codes[:, 0::4] = packed & 3
    codes[:, 1::4] = (packed >> 2) & 3
    codes[:, 2::4] = (packed >> 4) & 3
    codes[:, 3::4] = (packed >> 6) & 3
    return BED_CODE_TO_DOSAGE[codes[:, :n_samples]]


def write_bed(path: str | Path, codes: np.ndarray) -> None:
    """Write a .bed file from a (n_sites, n_samples) matrix of 2-bit codes."""
    n_sites, n_samples = codes.shape
    pad = (-n_samples) % 4
    if pad:
        filler = np.full((n_sites, pad), 1, dtype=np.uint8)  # code 01 = missing
        codes = np.concatenate([codes, filler], axis=1)
    groups = codes.astype(np.uint8).reshape(n_sites, -1, 4)
    packed = groups[:, :, 0] | (groups[:, :, 1] << 2) | (groups[:, :, 2] << 4) | (groups[:, :, 3] << 6)
    with open(path, "wb") as handle:
        handle.write(PLINK_MAGIC)
        handle.write(packed.tobytes())


def write_bim(path: str | Path, records: list[dict]) -> None:
    with open(path, "w") as handle:
        for record in records:
            handle.write(
                f"{record['chrom']}\t{record['rsid']}\t{record.get('cm', 0)}\t"
                f"{record['pos']}\t{record['a1']}\t{record['a2']}\n"
            )


def write_fam(path: str | Path, samples: list[dict]) -> None:
    with open(path, "w") as handle:
        for sample in samples:
            handle.write(
                f"{sample.get('fid', '0')}\t{sample['iid']}\t0\t0\t{sample.get('sex', '0')}\t-9\n"
            )
