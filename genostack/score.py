"""Significance filtering, goal weighting and grouping of findings into functional modules."""
from __future__ import annotations

from dataclasses import dataclass, field

from .annotate import Finding
from .config import EVIDENCE_WEIGHT
from .panel import Panel

MIN_FINAL_SCORE = 3.0          # auto findings below this are dropped
MINOR_PANEL_IMPACT = 1.0       # panel findings with impact 1 are kept as "minor" context
MAX_AUTO_PER_MODULE = 12
MAX_OTHER = 30

MODULE_ORDER_HINT = ["methylation", "neuro", "sleep_caffeine", "energy_mito", "vitamins_minerals", "inflammation_gut",
                     "cardiometabolic_longevity", "performance", "pharmaco", "other"]


@dataclass
class ModuleBlock:
    id: str
    title: str
    goals: list[str]
    description: str
    findings: list[Finding] = field(default_factory=list)   # major (sorted by score desc)
    minor: list[Finding] = field(default_factory=list)      # impact-1 panel context
    genes: list[str] = field(default_factory=list)
    weight: float = 0.0

    def all_findings(self) -> list[Finding]:
        return self.findings + self.minor


@dataclass
class Analysis:
    modules: list[ModuleBlock]
    stats: dict
    goals: list[str]
    n_raw_findings: int


def _goal_mult(module_goals: list[str], goals: list[str]) -> float:
    if not module_goals or not goals:
        return 1.0
    return 1.0 if set(module_goals) & set(goals) else 0.85


def rank(findings: list[Finding], panel: Panel, goals: list[str]) -> Analysis:
    gene_mod = panel.gene_module()
    mods: dict[str, ModuleBlock] = {
        mid: ModuleBlock(mid, meta["title"], list(meta.get("goals", [])), meta.get("description", ""))
        for mid, meta in panel.modules.items()
    }
    mods.setdefault("pharmaco", ModuleBlock("pharmaco", "Farmacogenómica", ["pharmaco", "health"], ""))
    mods["other"] = ModuleBlock("other", "Otros hallazgos significativos (cruce automático fuera del panel)", ["health"], "")

    panel_rsids = {f.rsid for f in findings if f.source in ("panel",)}
    covered_rsids = {v.rsid for v in panel.variants} | {r for rule in panel.rules for r in rule.rsids}  # judged by the panel
    haplo_rsids = {r for f in findings if f.source == "haplotype" for r in f.rsid.split(",")}
    # index auto findings by rsid to attach as support to panel findings
    support: dict[str, list[Finding]] = {}
    for f in findings:
        if f.source in ("snpedia", "clinvar", "gwas", "pgx") and f.rsid in panel_rsids:
            support.setdefault(f.rsid, []).append(f)

    for f in findings:
        if f.module is None:
            f.module = gene_mod.get(f.gene.upper(), "other") if f.gene else "other"
        ev = EVIDENCE_WEIGHT.get(f.evidence, 0.5)
        mult = _goal_mult(mods.get(f.module, mods["other"]).goals, goals)
        f.extra["final"] = round(f.score * (0.75 + 0.25 * ev) * mult, 2)

    for f in findings:
        blk = mods.get(f.module) or mods["other"]
        if f.source in ("panel", "haplotype"):
            if f.rsid in support:
                f.extra["support"] = [
                    {"source": s.source, "title": s.title, "detail": s.detail[:200], "score": s.score} for s in support[f.rsid]
                ]
            if f.impact >= 2 or f.source == "haplotype" and f.impact >= 1.5:
                blk.findings.append(f)
            elif f.impact >= MINOR_PANEL_IMPACT:
                blk.minor.append(f)
        else:
            if f.rsid in panel_rsids or (f.source != "pgx" and f.rsid in haplo_rsids):
                continue  # already represented by panel/haplotype finding
            if f.rsid in covered_rsids and not (f.source == "clinvar" and "atogénica" in f.title):
                continue  # panel judged this site as normal for this genotype; don't resurface it via GWAS/SNPedia
            if f.extra["final"] < MIN_FINAL_SCORE:
                continue
            blk.findings.append(f)

    out: list[ModuleBlock] = []
    for mid, blk in mods.items():
        blk.findings.sort(key=lambda x: -x.extra["final"])
        blk.minor.sort(key=lambda x: -x.extra["final"])
        # cap auto findings per module (keep all panel/haplotype ones)
        auto = [f for f in blk.findings if f.source not in ("panel", "haplotype")]
        cap = MAX_OTHER if mid == "other" else MAX_AUTO_PER_MODULE
        if len(auto) > cap:
            drop = set(id(f) for f in auto[cap:])
            blk.findings = [f for f in blk.findings if id(f) not in drop]
        blk.genes = list(dict.fromkeys(g for f in blk.all_findings() for g in [f.gene] if g))
        blk.weight = round(sum(f.extra["final"] for f in blk.findings) + 0.3 * sum(f.extra["final"] for f in blk.minor), 1)
        if blk.findings or blk.minor:
            out.append(blk)
    out.sort(key=lambda b: (-b.weight if b.id != "other" else 1e9, MODULE_ORDER_HINT.index(b.id) if b.id in MODULE_ORDER_HINT else 99))
    return Analysis(modules=out, stats={}, goals=goals, n_raw_findings=len(findings))
