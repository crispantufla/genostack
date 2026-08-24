"""Raw-report parsers: 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, minimal VCF.

Output: a Genome = dict rsid -> Call(chrom, pos, genotype) with genotype letters
sorted alphabetically (e.g. "AG"), no-calls dropped, alleles on GRCh37 plus strand
(which is what all consumer raw files provide).
"""
from __future__ import annotations

import csv
import gzip
import io
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path

VALID = set("ACGTDI")


@dataclass(frozen=True)
class Call:
    chrom: str
    pos: int
    genotype: str  # sorted letters, e.g. "AG", "A" (haploid), "DI"


@dataclass
class Genome:
    source_format: str
    build: str
    calls: dict[str, Call]
    n_total_lines: int
    n_nocalls: int

    @property
    def n_calls(self) -> int:
        return len(self.calls)

    def get(self, rsid: str) -> str | None:
        c = self.calls.get(rsid)
        return c.genotype if c else None

    def sex_guess(self) -> str:
        """Crude: many heterozygous X calls -> female; Y calls present -> male."""
        x_het = sum(1 for r, c in self.calls.items() if c.chrom == "X" and len(c.genotype) == 2 and c.genotype[0] != c.genotype[1])
        x_tot = sum(1 for c in self.calls.values() if c.chrom == "X")
        y_tot = sum(1 for c in self.calls.values() if c.chrom == "Y")
        if x_tot and x_het / x_tot > 0.05:
            return "female"
        if y_tot > 50:
            return "male"
        return "unknown"


def _open_text(path: Path) -> io.TextIOBase:
    data = path.read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    elif data[:2] == b"PK":
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            names = [n for n in z.namelist() if not n.endswith("/")]
            names.sort(key=lambda n: -z.getinfo(n).file_size)
            data = z.read(names[0])
    return io.StringIO(data.decode("utf-8", errors="replace"))


def _norm_chrom(c: str) -> str:
    c = c.strip().upper().removeprefix("CHR")
    return {"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"}.get(c, c)


def _norm_gt(gt: str) -> str | None:
    gt = gt.strip().upper().replace("/", "").replace("|", "")
    if not gt or gt in {"--", "00", "0", "NC", "-", "N"}:
        return None
    if any(ch not in VALID for ch in gt):
        return None
    return "".join(sorted(gt))


def parse_raw(path: str | Path) -> Genome:
    path = Path(path)
    fh = _open_text(path)
    text = fh.read()
    head = text[:4000].lower()
    if "23andme" in head or "# rsid\tchromosome\tposition\tgenotype" in head:
        return _parse_23andme(text)
    if "ancestrydna" in head or "allele1\tallele2" in head:
        return _parse_ancestry(text)
    if "##fileformat=vcf" in head:
        return _parse_vcf(text)
    if "myheritage" in head or head.startswith("rsid,chromosome,position,result") or '"rsid","chromosome","position","result"' in head:
        return _parse_myheritage(text)
    if "rsid,chromosome,position,result" in head or "rsid,chromosome,position,genotype" in head:
        return _parse_ftdna(text)
    # fallback: sniff by columns of first data line
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        if "\t" in line:
            parts = line.split("\t")
            if len(parts) >= 5:
                return _parse_ancestry(text)
            return _parse_23andme(text)
        if "," in line:
            return _parse_ftdna(text)
    raise ValueError("Formato de raw report no reconocido")


def _build_from_header(text: str) -> str:
    m = re.search(r"(GRCh3[78]|build\s*3[78]|hg1[89])", text[:5000], re.I)
    if not m:
        return "GRCh37"
    s = m.group(1).lower()
    return "GRCh38" if ("38" in s or "hg38" in s) else "GRCh37"


def _parse_23andme(text: str) -> Genome:
    calls: dict[str, Call] = {}
    n = nocall = 0
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\r").split("\t")
        if len(parts) < 4:
            continue
        rsid, chrom, pos, gt = parts[0], parts[1], parts[2], parts[3]
        n += 1
        g = _norm_gt(gt)
        if g is None:
            nocall += 1
            continue
        try:
            calls[rsid] = Call(_norm_chrom(chrom), int(pos), g)
        except ValueError:
            continue
    return Genome("23andMe", _build_from_header(text), calls, n, nocall)


def _parse_ancestry(text: str) -> Genome:
    calls: dict[str, Call] = {}
    n = nocall = 0
    for line in text.splitlines():
        if not line or line.startswith("#") or line.lower().startswith("rsid"):
            continue
        parts = line.rstrip("\r").split("\t")
        if len(parts) < 5:
            continue
        rsid, chrom, pos, a1, a2 = parts[:5]
        n += 1
        g = _norm_gt(a1 + a2)
        if g is None:
            nocall += 1
            continue
        try:
            calls[rsid] = Call(_norm_chrom(chrom), int(pos), g)
        except ValueError:
            continue
    return Genome("AncestryDNA", _build_from_header(text), calls, n, nocall)


def _parse_csvlike(text: str, fmt: str) -> Genome:
    calls: dict[str, Call] = {}
    n = nocall = 0
    reader = csv.reader(l for l in text.splitlines() if l and not l.startswith("#"))
    for row in reader:
        if len(row) < 4 or row[0].lower() in {"rsid", "rs_id"}:
            continue
        rsid, chrom, pos, gt = row[0].strip(), row[1], row[2], row[3]
        n += 1
        g = _norm_gt(gt)
        if g is None:
            nocall += 1
            continue
        try:
            calls[rsid] = Call(_norm_chrom(chrom), int(pos), g)
        except ValueError:
            continue
    return Genome(fmt, _build_from_header(text), calls, n, nocall)


def _parse_myheritage(text: str) -> Genome:
    return _parse_csvlike(text, "MyHeritage")


def _parse_ftdna(text: str) -> Genome:
    return _parse_csvlike(text, "FamilyTreeDNA")


def _parse_vcf(text: str) -> Genome:
    calls: dict[str, Call] = {}
    n = nocall = 0
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        chrom, pos, rsid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
        if not rsid.startswith("rs"):
            continue
        n += 1
        fmt = parts[8].split(":")
        sample = parts[9].split(":")
        try:
            gt_field = sample[fmt.index("GT")]
        except (ValueError, IndexError):
            nocall += 1
            continue
        alleles = [ref] + alt.split(",")
        idx = re.split(r"[/|]", gt_field)
        try:
            letters = [alleles[int(i)] for i in idx if i != "."]
        except (ValueError, IndexError):
            nocall += 1
            continue
        if not letters or any(len(a) != 1 for a in letters):
            nocall += 1
            continue
        calls[rsid] = Call(_norm_chrom(chrom), int(pos), "".join(sorted(letters)))
    return Genome("VCF", _build_from_header(text), calls, n, nocall)


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def complement(gt: str) -> str:
    return "".join(sorted(gt.translate(COMPLEMENT)))
