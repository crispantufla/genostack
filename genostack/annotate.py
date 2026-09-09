"""Cross the genome against the local DuckDB knowledge base and the curated panel.

Produces a flat list of `Finding`s (not yet filtered/ranked — see score.py).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import duckdb

from .config import (CLINVAR_MIN_STARS, DB_PATH, GWAS_MAX_PLAUSIBLE_OR, GWAS_MIN_BETA_SD, GWAS_MIN_OR,
                     GWAS_STRONG_OR, GWAS_STRONG_OR_MIN_N, GWAS_TRAIT_KEYWORDS, SNPEDIA_MIN_MAGNITUDE)
from .parsers import Genome, complement
from .panel import PanelVariant, HaploRule, load_panel

AMBIG = ({"A", "T"}, {"C", "G"})


@dataclass
class Finding:
    source: str            # panel | haplotype | snpedia | clinvar | gwas | pgx
    rsid: str              # may be a comma list for haplotypes
    gene: str
    genotype: str          # user's genotype (plus strand), or diplotype name for haplotypes
    title: str             # short headline
    detail: str            # phenotype / interpretation
    score: float           # 0-10 composite (pre goal-weighting)
    evidence: str          # A-E
    direction: str = "risk"  # risk | beneficial | mixed | info
    module: str | None = None
    impact: float = 0.0
    effect_size: str = ""
    mechanism: str = ""
    levers: list[str] = field(default_factory=list)
    caveats: str = ""
    refs: list[str] = field(default_factory=list)
    extra: dict = field(default_factory=dict)
    partial: bool = False

    def key(self) -> str:
        return f"{self.source}:{self.rsid}:{self.title}"


def _is_ambiguous(a: str, b: str) -> bool:
    return {a, b} in AMBIG


def _copies(genotype: str, allele: str) -> int:
    return sum(1 for ch in genotype if ch == allele)


def _sorted(g: str) -> str:
    return "".join(sorted(g))


def open_db() -> duckdb.DuckDBPyConnection:
    if not DB_PATH.exists():
        raise SystemExit(f"No existe la base local {DB_PATH}. Ejecuta primero: genostack fetch-data")
    return duckdb.connect(str(DB_PATH), read_only=True)


def load_genome_into(con: duckdb.DuckDBPyConnection, genome: Genome) -> None:
    import pyarrow as pa
    items = list(genome.calls.items())
    tbl = pa.table({"rsid": [r for r, _ in items], "genotype": [c.genotype for _, c in items], "chrom": [c.chrom for _, c in items]})
    con.register("user_gt_arrow", tbl)
    con.execute("CREATE OR REPLACE TEMP TABLE user_gt AS SELECT * FROM user_gt_arrow")
    con.unregister("user_gt_arrow")


# ----------------------------------------------------------------------------- panel
def eval_panel(genome: Genome, variants: list[PanelVariant]) -> list[Finding]:
    out: list[Finding] = []
    for v in variants:
        gt = genome.get(v.rsid)
        if gt is None:
            continue
        g = v.genotypes.get(gt)
        if g is None and len(gt) == 1:   # haploid call vs diploid keys
            g = v.genotypes.get(gt + gt)
        if g is None:
            continue
        impact = float(g.get("impact", 0))
        if impact <= 0:
            continue
        out.append(Finding(
            source="panel", rsid=v.rsid, gene=v.gene, genotype=gt,
            title=f"{v.gene} {v.name}".strip(), detail=g.get("phenotype", ""),
            score=min(10.0, impact * 2.2), evidence=v.evidence, direction=g.get("direction", v.direction),
            module=v.module, impact=impact, effect_size=v.effect_size, mechanism=v.mechanism,
            levers=list(v.levers), caveats=v.caveats, refs=list(v.refs),
            extra={"panel_id": v.id, "effect_allele": v.effect_allele, "copies": _copies(gt, v.effect_allele)},
        ))
    return out


def eval_haplotypes(genome: Genome, rules: list[HaploRule]) -> list[Finding]:
    out: list[Finding] = []
    for r in rules:
        res = r.evaluate(genome)
        if res is None:
            continue
        impact = float(res.get("impact", 0))
        if impact <= 0:
            continue
        out.append(Finding(
            source="haplotype", rsid=",".join(r.rsids), gene=r.gene, genotype=res.get("name", ""),
            title=f"{r.title}: {res.get('name','')}",
            detail=res.get("phenotype") or res.get("description") or (r.mechanism[:260] if r.mechanism else ""),
            score=min(10.0, impact * 2.2), evidence=r.evidence, direction=res.get("direction", "risk"),
            module=r.module, impact=impact, effect_size=res.get("effect_size", ""), mechanism=r.mechanism,
            levers=list(res.get("levers") or r.levers), caveats=r.caveats, refs=list(r.refs),
            extra={"rule_id": r.id, "genotypes": res.get("genotypes", {})}, partial=bool(res.get("partial")),
        ))
    return out


# ----------------------------------------------------------------------------- SNPedia
def eval_snpedia(con: duckdb.DuckDBPyConnection) -> list[Finding]:
    rows = con.execute("""
        SELECT g.rsid, u.genotype, g.allele1, g.allele2, g.magnitude, g.repute, g.summary,
               s.gene, s.orientation, s.summary, u.chrom
        FROM user_gt u
        JOIN snpedia_geno g ON g.rsid = u.rsid
        LEFT JOIN snpedia_snp s ON s.rsid = u.rsid
        WHERE g.magnitude >= ?
    """, [SNPEDIA_MIN_MAGNITUDE]).fetchall()
    out: list[Finding] = []
    seen: set[str] = set()
    for rsid, ugt, a1, a2, mag, repute, gsum, gene, orient, ssum, chrom in rows:
        a1, a2 = (a1 or "").upper(), (a2 or "").upper()
        if not a1 or not a2 or len(a1) != 1 or len(a2) != 1:
            continue
        sg = _sorted(a1 + a2)
        if orient == "minus":
            sg = complement(sg)
        ug = ugt if len(ugt) == 2 else ugt + ugt
        if sg != ug:
            continue
        if rsid in seen:
            continue
        seen.add(rsid)
        direction = {"good": "beneficial", "bad": "risk"}.get(repute or "", "mixed")
        ev = "B" if mag >= 4 else "C"
        out.append(Finding(
            source="snpedia", rsid=rsid, gene=gene or "", genotype=ugt,
            title=f"{gene or rsid} — SNPedia magnitud {mag:g}" if gene else f"{rsid} — SNPedia magnitud {mag:g}",
            detail=(gsum or ssum or "").strip(), score=min(10.0, float(mag)), evidence=ev, direction=direction,
            impact=min(4.0, float(mag) / 2.5), extra={"magnitude": mag, "repute": repute, "snp_summary": ssum},
            caveats="Orientación SNPedia " + (orient or "desconocida") + "; magnitud es una escala subjetiva de SNPedia (0-10).",
        ))
    return out


# ----------------------------------------------------------------------------- ClinVar
def eval_clinvar(con: duckdb.DuckDBPyConnection) -> list[Finding]:
    rows = con.execute("""
        SELECT c.rsid, u.genotype, c.gene, c.name, c.clinsig, c.stars, c.review_status, c.phenotypes, c.ref, c.alt, c.vtype, c.variation_id, c.n_submitters
        FROM user_gt u JOIN clinvar c ON c.rsid = u.rsid
        WHERE c.stars >= ?
    """, [CLINVAR_MIN_STARS]).fetchall()
    out: list[Finding] = []
    best: dict[str, Finding] = {}
    n_indel_skipped = 0
    for rsid, ugt, gene, name, clinsig, stars, review, phen, ref, alt, vtype, vid, nsub in rows:
        alt = (alt or "").upper()
        ref = (ref or "").upper()
        if set(ugt) <= {"D", "I"}:
            # Consumer arrays code indel sites as D/I relative to the probe's own definition, which says
            # nothing about WHICH indel is present: several distinct ClinVar indels share one rsID, and the
            # longer allele ("I") is often the reference. Matching them produced homozygous "pathogenic"
            # calls in CDKL5/MECP2/F8 for a healthy adult — exactly backwards. They are not interpretable.
            n_indel_skipped += 1
            continue
        if vtype and vtype not in ("single nucleotide variant",):
            continue          # an indel/CNV/microsatellite cannot be read off an A/C/G/T array call
        if len(alt) != 1 or len(ref) != 1:
            continue
        copies = _copies(ugt, alt)
        if copies == 0:
            continue
        cs = (clinsig or "").lower()
        if "pathogenic" in cs and "conflicting" not in cs:
            base = 6.0 if copies == 1 else 8.0
            direction = "risk"
            kind = "Patogénica" if "likely" not in cs else "Probablemente patogénica"
        elif "risk factor" in cs:
            base, direction, kind = 4.5, "risk", "Factor de riesgo"
        elif "drug response" in cs:
            base, direction, kind = 5.0, "info", "Respuesta a fármaco"
        elif "protective" in cs:
            base, direction, kind = 4.5, "beneficial", "Protectora"
        elif "affects" in cs:
            base, direction, kind = 4.0, "mixed", "Afecta"
        else:
            continue
        score = min(10.0, base + 0.5 * (stars - 2))
        zyg = "homocigoto/hemicigoto" if copies >= 2 or len(ugt) == 1 else "heterocigoto (portador)"
        amb = _is_ambiguous(ref, alt) if len(alt) == 1 and len(ref) == 1 else False
        f = Finding(
            source="clinvar", rsid=rsid, gene=gene or "", genotype=ugt,
            title=f"{gene or rsid}: {kind} ({zyg})", detail=f"{name} — {clinsig}; {phen}"[:400],
            score=score, evidence="A" if stars >= 3 else "B", direction=direction, impact=min(4.0, score / 2.2),
            extra={"stars": stars, "review": review, "variation_id": vid, "copies": copies, "clinsig": clinsig},
            caveats="Los arrays de consumo tienen alta tasa de falsos positivos en variantes raras: confirmar con secuenciación clínica antes de actuar."
                    + (" SNP A/T o C/G: orientación de hebra no resoluble con certeza." if amb else ""),
        )
        if rsid not in best or f.score > best[rsid].score:
            best[rsid] = f
    eval_clinvar.n_indel_skipped = n_indel_skipped   # surfaced in the report's limitations
    return list(best.values())


# ----------------------------------------------------------------------------- GWAS
_KW_RE = re.compile("|".join(re.escape(k) for k in GWAS_TRAIT_KEYWORDS), re.I)
_BLOCK_RE = re.compile(r"protein level|protein measurement|gene expression|splic|methylation|particle|very large|very small|"
                       r"\bratio\b|metabolite|concentration of|in medium|in large|in small|lipoprotein subclass|"
                       r"cholesteryl esters|free cholesterol|phospholipids in|triglycerides in|apolipoprotein a1 in|"
                       r"esterified|diameter|imaging|brain morphology|cortical|white matter|electrocardio|"
                       r"\bmean\b.*\bvolume\b|intensity|\bscan\b|phosphatidyl|sphingomyelin|ceramide|lysophosph|acylcarnitine|"
                       r"glycerophospho|\(\d+:\d+|level in |levels in |interaction|df test|\bx\b|xenobiotic|adjusted for|"
                       r"conditional|pleiotrop|multivariate|joint analysis|\bPC\(|\bTG\(|\bLPC\(|\bSM\(|\bCE\(|"
                       r"plasma protein|serum protein|levels of protein|protein levels|antigen|urinary|\bpeptide\b|response to|"
                       r"treatment|survival|change in|trajectory|longitudinal|\bin .* patients\b|\bin .* cases\b|cases vs|\bamong\b|"
                       r"triacylglycerol|diacylglycerol|to total lipids|\(model|\bfa\d|\bpc aa|\bpc ae|lysopc|\bsm c\d|"
                       r"time to event|metastasis|mortality|prognosis|progression|recurrence|\bstage\b|\bin (type [12] )?diabet", re.I)
# traits whose level is directly actionable via supplementation/lifestyle even with small per-allele betas
_ACTIONABLE_RE = re.compile(r"vitamin|folate|b12|cobalamin|homocysteine|ferritin|\biron\b|transferrin|hemoglobin|haemoglobin|zinc|"
                            r"magnesium|selenium|copper|urate|uric|omega|docosahexaenoic|eicosapentaenoic|arachidonic|linoleic|"
                            r"choline|carnitine|glutathione|25-hydroxy|retinol|carotene|tocopherol|testosterone|shbg|\bigf-?1\b|"
                            r"c-reactive|\bcrp\b|\btsh\b|thyro|caffeine|coffee|sleep|chronotype|hba1c|glucose|insulin|triglyceride|"
                            r"\bldl\b|\bhdl\b|cholesterol|apolipoprotein [ab]\b|lipoprotein\(a\)|lp\(a\)|blood pressure|bone mineral|"
                            r"creatinine|egfr|cystatin|alanine aminotransferase|gamma-glutamyl|bilirubin|estradiol|cortisol|"
                            r"dhea|melatonin|ferritin|hepcidin|b6|pyridox|riboflavin|thiamin|niacin|biotin|iodine|albumin\b", re.I)


def _sample_n(text: str | None) -> int:
    if not text:
        return 0
    nums = [int(x.replace(",", "")) for x in re.findall(r"\d[\d,]{2,}", text)]
    return sum(nums) if nums else 0


def _beta_keep(ob: float, ci: str, trait: str) -> tuple[bool, float]:
    """Return (keep, magnitude_score) for beta-type associations, unit-aware."""
    ci_l = (ci or "").lower()
    ab = abs(ob)
    if re.search(r"\bsd\b|z-score|z score|standard deviation|s\.d\.|\blog\b|\bln\b", ci_l):
        if ab < GWAS_MIN_BETA_SD or ab > 1.5:   # >1.5 SD per allele is implausible (catalog artefacts)
            return False, 0.0
        return True, min(10.0, 3.0 + 10.0 * min(ab, 0.7))
    if "%" in ci_l:
        if ab < 5 or ab > 60:
            return False, 0.0
        return True, min(10.0, 3.0 + ab / 5.0)
    # physical/unknown units (mg/dL, cm, kg, years, "unit"...) — unknown scale
    if _ACTIONABLE_RE.search(trait):
        return True, 3.2  # actionable biomarker: keep, modest score (ranking then via p/N)
    return False, 0.0


def eval_gwas(con: duckdb.DuckDBPyConnection) -> list[Finding]:
    # 47 % of catalog rows carry no risk-allele frequency, which used to blind the major-allele filter below
    # (a "risk" allele present in 92 % of people was reported as a finding). Fill it in from the other rows
    # of the same rsID + allele, which is the same quantity measured by a study that did report it.
    rows = con.execute("""
        WITH freq AS (
            SELECT rsid, risk_allele, median(raf) AS raf_med
            FROM gwas WHERE raf IS NOT NULL GROUP BY 1, 2
        )
        SELECT g.rsid, u.genotype, g.risk_allele, g.trait, g.mapped_trait, g.pvalue, g.or_beta, g.effect_type,
               g.ci_text, g.gene, g.pubmedid, g.first_author, g.pub_date, g.sample_text,
               COALESCE(g.raf, f.raf_med) AS raf
        FROM user_gt u
        JOIN gwas g ON g.rsid = u.rsid
        LEFT JOIN freq f ON f.rsid = g.rsid AND f.risk_allele = g.risk_allele
        WHERE g.or_beta IS NOT NULL AND g.risk_allele IN ('A','C','G','T')
    """).fetchall()
    best: dict[tuple[str, str], tuple[tuple, Finding]] = {}
    for rsid, ugt, ra, trait, mtrait, p, ob, etype, ci, gene, pmid, author, date, stext, raf in rows:
        tname = (mtrait or trait or "")
        if not (_KW_RE.search(tname) or _KW_RE.search(trait or "")):
            continue
        if _BLOCK_RE.search(trait or "") or _BLOCK_RE.search(tname):
            continue
        if re.search(r"\blevels?\b|measurement", trait or "", re.I) and not _ACTIONABLE_RE.search(trait or ""):
            continue  # molecular QTL traits (proteins, lipid species...) unless an actionable biomarker
        if etype == "or":
            if ob <= 0:
                continue
            eff = 1 / ob if ob < 1 else ob
            if eff < GWAS_MIN_OR:
                continue
            # No common variant moves a complex trait this much. Values like a flat OR=100 repeated across
            # dozens of rsIDs of one study are data-entry artefacts, and they headline the report if kept.
            if eff > GWAS_MAX_PLAUSIBLE_OR:
                continue
            if eff > GWAS_STRONG_OR and _sample_n(stext) < GWAS_STRONG_OR_MIN_N:
                continue   # a genuinely strong effect needs a well-powered study behind it
            mag_score = min(10.0, 3.0 + 2.0 * math.log2(eff))
            eff_text = f"OR {ob:g} por alelo {ra}" + (" (protector)" if ob < 1 else "")
        else:
            keep, mag_score = _beta_keep(ob, ci, f"{trait} {tname}")
            if not keep:
                continue
            eff_text = f"β {ob:g} {ci or ''}".strip()
        copies = _copies(ugt, ra)
        if copies == 0:
            continue
        if raf is not None and raf >= 0.5:
            continue  # effect allele is the major allele: carrying it is the population-typical state, not a finding
        amb = {ra, complement(ra)} in AMBIG and set(ugt) <= {ra, complement(ra)}
        n = _sample_n(stext)
        if n < 5000 and not (etype == "or" and ob >= 2.0):
            continue  # small studies: only keep very large ORs
        bonus = min(1.0, (math.log10(1 / p) if p and p > 0 else 8) / 100.0 + (0.3 if n > 100_000 else 0.0))
        rank = (n, -float(p or 1))
        f = Finding(
            source="gwas", rsid=rsid, gene=(gene or "").split(" - ")[0].split(",")[0].strip(), genotype=ugt,
            title=f"{trait}", detail=f"{eff_text}; p={p:.1e}; {copies} copia(s) del alelo {ra}; {author} {date[:4] if date else ''} (PMID {pmid}); N≈{n:,}",
            score=min(10.0, (mag_score + bonus) * (1.15 if copies == 2 else 1.0)), evidence="B",
            direction=("beneficial" if (etype == "or" and ob < 1) else "risk" if etype == "or" or ob > 0 else "mixed"), impact=min(4.0, mag_score / 2.2),
            extra={"risk_allele": ra, "copies": copies, "pvalue": p, "or_beta": ob, "effect_type": etype, "pmid": pmid, "n": n, "raf": raf, "mapped_trait": mtrait},
            caveats=("SNP A/T o C/G: orientación del alelo de riesgo no verificable con certeza. " if amb else "") +
                    "Asociación GWAS (efecto por alelo a nivel poblacional); no es un diagnóstico.",
        )
        k = (rsid, (mtrait or trait).lower())
        if k not in best or rank > best[k][0]:
            best[k] = (rank, f)
    # cap traits per rsid (pleiotropic loci otherwise flood the report)
    per_rsid: dict[str, list[Finding]] = {}
    for _, f in best.values():
        per_rsid.setdefault(f.rsid, []).append(f)
    out: list[Finding] = []
    for fs in per_rsid.values():
        fs.sort(key=lambda x: (-x.score, -x.extra.get("n", 0)))
        out.extend(fs[:3])
    return out


# ----------------------------------------------------------------------------- Pharmacogenomics (Open Targets / PharmGKB)
# Many PharmGKB annotations describe the *reference* state ("do not carry a copy of...", "assigned normal
# function"). Reporting them is noise: they say the person is normal, and they crowded out the real hits.
_PGX_NORMAL_RE = re.compile(r"do(?:es)? not (?:have|carry)|assigned normal function|"
                            r"normal (?:function|metaboli[sz]|activity) (?:by|as compared)", re.I)
def _pgx_ref_allele(comparisons: list[str] | None) -> str | None:
    """Most frequent single-letter 'comparison' allele across annotations ~ reference allele."""
    from collections import Counter
    c: Counter = Counter()
    for x in (comparisons or []):
        x = (x or "").upper().strip()
        if len(x) == 1 and x in "ACGT":
            c[x] += 2
        elif len(x) == 2 and x[0] == x[1] and x[0] in "ACGT":
            c[x[0]] += 1
    return c.most_common(1)[0][0] if c else None


def eval_pgx(con: duckdb.DuckDBPyConnection) -> list[Finding]:
    rows = con.execute("""
        SELECT p.rsid, u.genotype, p.genotype, p.annotation, p.level, p.category, p.phenotype, p.drugs, p.literature, t.symbol,
               p.comparisons
        FROM user_gt u JOIN ot_pgx p ON p.rsid = u.rsid
        LEFT JOIN ot_target t ON t.ensg = p.ensg
        WHERE p.level IN ('1A','1B','2A','2B') AND p.genotype IS NOT NULL
    """).fetchall()
    comp_by_rsid: dict[str, list] = {}
    for r in rows:
        comp_by_rsid.setdefault(r[0], []).extend(r[10] or [])
    ref = {rs: _pgx_ref_allele(v) for rs, v in comp_by_rsid.items()}
    grouped: dict[str, dict] = {}
    for rsid, ugt, pgt, ann, level, cat, phen, drugs, lit, sym, comps in rows:
        pg = re.sub(r"[^A-Z]", "", (pgt or "").upper())
        if len(pg) not in (1, 2):
            continue
        ug = ugt if len(ugt) == 2 else ugt + ugt
        if _sorted(pg if len(pg) == 2 else pg + pg) != ug:
            continue
        ra = ref.get(rsid)
        if ra and all(ch == ra for ch in ug):
            continue  # homozygous reference: the annotation describes the normal state
        if _PGX_NORMAL_RE.search(ann or phen or ""):
            continue  # the text itself says the person lacks the variant / has normal function
        g = grouped.setdefault(rsid, {"gene": sym or "", "gt": ugt, "levels": set(), "drugs": [], "texts": [], "cats": set(), "pmids": set()})
        g["levels"].add(level)
        for d in (drugs or []):
            if d not in g["drugs"]:
                g["drugs"].append(d)
        g["cats"].add(cat or "")
        t = (ann or phen or "").strip()
        if t and t not in g["texts"]:
            g["texts"].append(t)
        for x in (lit or []):
            g["pmids"].add(x)
    out: list[Finding] = []
    for rsid, g in grouped.items():
        best_level = sorted(g["levels"])[0]
        score = {"1A": 7.5, "1B": 7.0, "2A": 5.5, "2B": 5.0}[best_level]
        out.append(Finding(
            source="pgx", rsid=rsid, gene=g["gene"], genotype=g["gt"],
            title=f"{g['gene'] or rsid} × {', '.join(g['drugs'][:6])}{'…' if len(g['drugs']) > 6 else ''}",
            detail=" | ".join(g["texts"][:3])[:600],
            score=score, evidence="A" if best_level.startswith("1") else "B", direction="info", impact=min(4.0, score / 2.2),
            module="pharmaco",
            extra={"level": best_level, "categories": sorted(g["cats"]), "drugs": g["drugs"], "pmids": sorted(g["pmids"])[:10]},
            caveats="Anotación clínica PharmGKB/ClinPGx nivel " + best_level + ".",
        ))
    return out


# ----------------------------------------------------------------------------- gene context for the agent
def gene_context(con: duckdb.DuckDBPyConnection, genes: list[str], max_drugs: int = 25) -> dict[str, dict]:
    """Function text, known/experimental drugs (Open Targets, DGIdb) and CPIC pairs per gene."""
    genes = [g for g in dict.fromkeys(genes) if g]
    if not genes:
        return {}
    ctx: dict[str, dict] = {g: {} for g in genes}
    ph = ",".join("?" * len(genes))
    for sym, ensg, fn, pathways, tract in con.execute(
            f"SELECT symbol, ensg, function_text, pathways, tractability FROM ot_target WHERE symbol IN ({ph})", genes).fetchall():
        ctx[sym].update({"ensg": ensg, "function": (fn or "")[:700], "pathways": (pathways or [])[:8], "tractability": (tract or [])[:6]})
    ensg_to_sym = {v["ensg"]: k for k, v in ctx.items() if v.get("ensg")}
    if ensg_to_sym:
        eph = ",".join("?" * len(ensg_to_sym))
        rows = con.execute(f"""
            SELECT ct.ensg, m.name, m.drug_type, ct.max_stage, mo.action_type, mo.mechanism, ct.diseases
            FROM ot_clinical_target ct JOIN ot_molecule m ON m.chembl_id = ct.chembl_id
            LEFT JOIN ot_moa mo ON mo.chembl_id = ct.chembl_id AND list_contains(mo.target_ids, ct.ensg)
            WHERE ct.ensg IN ({eph})
            ORDER BY ct.ensg, ct.max_stage DESC""", list(ensg_to_sym)).fetchall()
        for ensg, name, dtype, stage, action, mech, dis in rows:
            sym = ensg_to_sym[ensg]
            lst = ctx[sym].setdefault("ot_drugs", [])
            if len(lst) < max_drugs:
                lst.append({"drug": name, "type": dtype, "max_stage": stage, "action": action, "mechanism": mech, "indications": (dis or [])[:4]})
        rows = con.execute(f"""
            SELECT mo.target_ids, m.name, m.max_stage, mo.action_type, mo.mechanism
            FROM ot_moa mo JOIN ot_molecule m ON m.chembl_id = mo.chembl_id
            WHERE list_has_any(mo.target_ids, ?)""", [list(ensg_to_sym)]).fetchall()
        for tids, name, stage, action, mech in rows:
            for t in tids:
                if t in ensg_to_sym:
                    lst = ctx[ensg_to_sym[t]].setdefault("ot_moa", [])
                    if len(lst) < max_drugs and not any(x["drug"] == name for x in lst):
                        lst.append({"drug": name, "max_stage": stage, "action": action, "mechanism": mech})
    for gene, drug, itype, score, approved in con.execute(
            f"SELECT gene, drug, interaction_type, max(score), bool_or(approved) FROM dgidb WHERE gene IN ({ph}) GROUP BY 1,2,3 ORDER BY 4 DESC NULLS LAST", genes).fetchall():
        lst = ctx[gene].setdefault("dgidb", [])
        if len(lst) < max_drugs:
            lst.append({"drug": drug, "type": itype, "approved": approved})
    for gene, drug, lvl, testing in con.execute(
            f"SELECT gene, drug, cpic_level, pgx_testing FROM cpic_pair WHERE gene IN ({ph}) AND cpic_level IN ('A','B')", genes).fetchall():
        ctx[gene].setdefault("cpic", []).append({"drug": drug, "level": lvl, "testing": testing})
    return ctx


# ----------------------------------------------------------------------------- orchestrate
def resolve_aliases(genome: Genome, con: duckdb.DuckDBPyConnection) -> int:
    """Map 23andMe internal ids (i3002432 ...) to rsIDs using SNPedia's alias table; returns number added."""
    has = con.execute("SELECT count(*) FROM information_schema.tables WHERE table_name='snpedia_alias'").fetchone()[0]
    if not has:
        return 0
    iids = [r for r in genome.calls if r.startswith("i") and r[1:].isdigit()]
    if not iids:
        return 0
    con.execute("CREATE OR REPLACE TEMP TABLE user_iids(iid VARCHAR)")
    con.executemany("INSERT INTO user_iids VALUES (?)", [(i,) for i in iids])
    n = 0
    for iid, rsid in con.execute("SELECT a.iid, a.rsid FROM snpedia_alias a JOIN user_iids u ON u.iid = a.iid").fetchall():
        if rsid not in genome.calls:
            genome.calls[rsid] = genome.calls[iid]
            n += 1
    return n


def annotate(genome: Genome, con: duckdb.DuckDBPyConnection) -> tuple[list[Finding], dict]:
    panel = load_panel()
    n_alias = resolve_aliases(genome, con)
    load_genome_into(con, genome)
    findings: list[Finding] = []
    findings += eval_panel(genome, panel.variants)
    findings += eval_haplotypes(genome, panel.rules)
    has_snpedia = con.execute("SELECT count(*) FROM information_schema.tables WHERE table_name='snpedia_geno'").fetchone()[0] > 0
    if has_snpedia:
        findings += eval_snpedia(con)
    findings += eval_clinvar(con)
    findings += eval_gwas(con)
    findings += eval_pgx(con)
    stats = {
        "n_calls": genome.n_calls,
        "n_in_snpedia": con.execute("SELECT count(DISTINCT u.rsid) FROM user_gt u JOIN snpedia_snp s ON s.rsid=u.rsid").fetchone()[0] if has_snpedia else 0,
        "n_in_clinvar": con.execute("SELECT count(DISTINCT u.rsid) FROM user_gt u JOIN clinvar c ON c.rsid=u.rsid").fetchone()[0],
        "n_in_gwas": con.execute("SELECT count(DISTINCT u.rsid) FROM user_gt u JOIN gwas g ON g.rsid=u.rsid").fetchone()[0],
        "n_in_pgx": con.execute("SELECT count(DISTINCT u.rsid) FROM user_gt u JOIN ot_pgx p ON p.rsid=u.rsid").fetchone()[0],
        "panel_variants": len(panel.variants), "panel_rules": len(panel.rules),
        "panel_covered": sum(1 for v in panel.variants if genome.get(v.rsid) is not None),
        "snpedia_available": has_snpedia,
        "n_alias_resolved": n_alias,
        "n_clinvar_indels_skipped": getattr(eval_clinvar, "n_indel_skipped", 0),
    }
    return findings, stats
