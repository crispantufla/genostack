"""Curated panel loader: variants + haplotype/combination rules (see SCHEMA.md)."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

PANEL_DIR = Path(__file__).resolve().parent


def _norm_gt(g: str) -> str:
    return "".join(sorted(str(g).strip().upper()))


@dataclass
class PanelVariant:
    id: str
    module: str
    gene: str
    rsid: str
    name: str
    effect_allele: str
    other_allele: str
    direction: str
    genotypes: dict[str, dict]
    effect_size: str
    evidence: str
    mechanism: str
    levers: list[str]
    caveats: str = ""
    refs: list[str] = field(default_factory=list)


@dataclass
class HaploRule:
    id: str
    module: str
    gene: str
    title: str
    type: str                       # combo | allele_count
    rsids: list[str]
    combos: dict[str, dict]         # for combo
    alleles: list[dict]             # for allele_count
    phenotypes: list[dict]          # for allele_count
    evidence: str
    mechanism: str
    levers: list[str]
    caveats: str = ""
    refs: list[str] = field(default_factory=list)

    def evaluate(self, genome) -> dict | None:
        gts = {r: genome.get(r) for r in self.rsids}
        present = {r: g for r, g in gts.items() if g}
        if not present:
            return None
        partial = len(present) < len(self.rsids)
        if self.type == "combo":
            if partial:
                return None  # combos need all sites
            for key, res in self.combos.items():
                conds = dict(p.split("=") for p in key.replace(" ", "").split(";"))
                if all(_norm_gt(v) == _norm_gt(gts[r]) or (len(gts[r]) == 1 and _norm_gt(v) == gts[r] * 2)
                       for r, v in conds.items() if r in gts):
                    out = dict(res)
                    out["genotypes"] = dict(present)
                    return out
            return None
        # allele_count
        loss = gain = count = 0
        for a in self.alleles:
            g = gts.get(a["rsid"])
            if not g:
                continue
            n = sum(1 for ch in g if ch == str(a["allele"]).upper())
            if len(g) == 1:
                n = 2 * n  # haploid call counts as homozygous
            fn = str(a.get("function", "loss")).lower()
            if fn in ("loss", "decreased", "slow", "risk", "reduced", "null"):
                loss += n
            elif fn in ("gain", "increased", "fast", "protective"):
                gain += n
            count += n
        env = {"loss": loss, "gain": gain, "count": count, "true": True}
        for ph in self.phenotypes:
            cond = str(ph.get("when", "true")).strip()
            if _eval_cond(cond, env):
                out = dict(ph)
                out["genotypes"] = dict(present)
                out["partial"] = partial
                out["detail_counts"] = {"loss": loss, "gain": gain, "count": count}
                return out
        return None


_COND_RE = re.compile(r"^(loss|gain|count)\s*(==|>=|<=|>|<)\s*(\d+)$")


def _eval_cond(cond: str, env: dict) -> bool:
    cond = cond.strip()
    if cond in ("true", "True", "1", ""):
        return True
    parts = [p.strip() for p in re.split(r"\band\b|&&", cond)]
    for p in parts:
        m = _COND_RE.match(p)
        if not m:
            return False
        v = env[m.group(1)]
        n = int(m.group(3))
        op = m.group(2)
        ok = {"==": v == n, ">=": v >= n, "<=": v <= n, ">": v > n, "<": v < n}[op]
        if not ok:
            return False
    return True


@dataclass
class Panel:
    variants: list[PanelVariant]
    rules: list[HaploRule]
    modules: dict[str, dict]   # module id -> {title, goals, description}

    def gene_module(self) -> dict[str, str]:
        m: dict[str, str] = {}
        for v in self.variants:
            m.setdefault(v.gene.upper(), v.module)
        for r in self.rules:
            m.setdefault(r.gene.upper(), r.module)
        return m


def load_panel(panel_dir: Path = PANEL_DIR) -> Panel:
    variants: list[PanelVariant] = []
    rules: list[HaploRule] = []
    modules: dict[str, dict] = {}
    seen_ids: set[str] = set()
    for f in sorted(panel_dir.glob("*.yaml")):
        data = yaml.safe_load(f.read_text(encoding="utf-8")) or {}
        if f.name == "haplotypes.yaml":
            for r in data.get("rules", []):
                rules.append(HaploRule(
                    id=r["id"], module=r.get("module") or _guess_module(r["gene"]), gene=r["gene"], title=r.get("title", r["id"]),
                    type=r["type"], rsids=[str(x) for x in r["rsids"]] if "rsids" in r else [a["rsid"] for a in r.get("alleles", [])],
                    combos=r.get("combos", {}), alleles=r.get("alleles", []), phenotypes=r.get("phenotypes", []),
                    evidence=str(r.get("evidence", "B")), mechanism=r.get("mechanism", ""), levers=list(r.get("levers", [])),
                    caveats=r.get("caveats", ""), refs=[str(x) for x in r.get("refs", [])]))
            continue
        mod = data.get("module") or f.stem
        modules[mod] = {"title": data.get("title", mod), "goals": data.get("goals", []), "description": data.get("description", "")}
        for v in data.get("variants", []):
            if v["id"] in seen_ids:
                raise ValueError(f"Panel: id duplicado {v['id']} en {f.name}")
            seen_ids.add(v["id"])
            gts = {_norm_gt(k): dict(val) for k, val in (v.get("genotypes") or {}).items()}
            variants.append(PanelVariant(
                id=v["id"], module=mod, gene=v["gene"], rsid=str(v["rsid"]).lower(), name=v.get("name", ""),
                effect_allele=str(v["effect_allele"]).upper(), other_allele=str(v.get("other_allele", "")).upper(),
                direction=v.get("direction", "risk"), genotypes=gts, effect_size=v.get("effect_size", ""),
                evidence=str(v.get("evidence", "C")), mechanism=v.get("mechanism", ""), levers=list(v.get("levers", [])),
                caveats=v.get("caveats", "") or "", refs=[str(x) for x in v.get("refs", [])]))
    # assign modules to rules lacking one: by gene, else by the module of any of their rsids
    gm = {v.gene.upper(): v.module for v in variants}
    rm = {v.rsid: v.module for v in variants}
    for r in rules:
        if r.module == "?":
            genes = [g.strip().upper() for g in re.split(r"[/,\s]+", r.gene) if g.strip()]
            mod = next((gm[g] for g in genes if g in gm), None) or next((rm[x] for x in r.rsids if x in rm), None)
            r.module = mod or "other"
    return Panel(variants, rules, modules)


_GUESS = {
    "APOE": "cardiometabolic_longevity", "CYP2C19": "pharmaco", "CYP2C9": "pharmaco", "CYP2D6": "pharmaco",
    "SLCO1B1": "pharmaco", "DPYD": "pharmaco", "TPMT": "pharmaco", "NAT2": "pharmaco", "MTHFR": "methylation",
    "HFE": "vitamins_minerals", "VDR": "vitamins_minerals", "GC": "vitamins_minerals", "FADS1": "vitamins_minerals",
    "FADS2": "vitamins_minerals", "COMT": "neuro", "ACTN3": "performance", "ACE": "performance", "VKORC1": "pharmaco",
    "UGT1A1": "pharmaco", "G6PD": "pharmaco", "CYP1A2": "sleep_caffeine",
}


def _guess_module(gene: str) -> str:
    return _GUESS.get(gene.upper(), "?")
