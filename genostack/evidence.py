"""Deterministic evidence audit for compounds proposed by the agent (Europe PMC + ClinicalTrials.gov).

No LLM involved: counts of human RCTs / meta-analyses, recent (last 5 y) human studies, and active trials,
plus a few top titles. Used to anchor/grade the agent's proposals.
"""
from __future__ import annotations

import asyncio
import re
from dataclasses import dataclass, field
from datetime import date

import httpx

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
CTG = "https://clinicaltrials.gov/api/v2/studies"


@dataclass
class EvidenceAudit:
    compound: str
    term: str = ""                   # actual search term used (the compound name is often a phrase)
    rct_meta_count: int = 0          # RCTs + meta-analyses + systematic reviews (any year)
    human_recent_count: int = 0      # human studies last 5 years
    total_count: int = 0             # any publication mentioning it in title
    active_trials: int = 0
    top_titles: list[str] = field(default_factory=list)
    trial_titles: list[str] = field(default_factory=list)
    error: str | None = None

    def grade_hint(self) -> str:
        if self.error:
            return "?"
        if self.rct_meta_count >= 20:
            return "A"
        if self.rct_meta_count >= 5:
            return "B"
        if self.human_recent_count >= 5 or self.rct_meta_count >= 1:
            return "C"
        if self.total_count >= 3:
            return "D"
        return "E"


_QUALIFIER_RE = re.compile(
    r"\b(?:very\s+)?(?:low|high|standard|micro|mega|full|half)[- ]?dose[sd]?\b|\bdose[sd]?\b|"
    r"supplementation|supplement|protocol|therap(?:y|ies)|treatments?|agents?|surveillance|avoidance|repletion|"
    r"preformed|re-?esterified|standardi[sz]ed|extract|monohydrate|bitartrate|glycinate|bisglycinate|"
    r"malate|citrate|chelate|liposomal|micellar|sublingual|topical|oral|intravenous|subcutaneous|"
    r"daily|nightly|weekly|as needed|targeted|adjunct|add-?on|combination|regimen|strategy|"
    r"instead of .*|in place of .*|rather than .*|for .*|when .*|if .*|plan\b|approach\b|"
    r"first[- ]line|second[- ]line|rescue|trial of|cycle[sd]?", re.I)
_SPLIT_RE = re.compile(r"\s*(?:/|\+|,|;|\band\b|\bor\b|\bplus\b|\bwith\b|\bvs\.?\b|—|–)\s*", re.I)
_GENERIC = {"iron", "blood", "donation", "food", "diet", "exercise", "sleep", "protein", "fat", "fiber",
            "fibre", "water", "salt", "sugar", "light", "training", "acid", "oil", "vitamin", "mineral",
            "day", "week", "night", "morning", "meal", "intake", "level", "levels", "test", "testing"}
# Drug-class / modality words. A candidate made ONLY of these is a modality, not a compound: searching it
# ("siRNA", "peptide therapy") returns the whole literature of the class and would grade an experimental
# molecule as if it had decades of trials behind it.
_CLASS = {"sirna", "rnai", "antisense", "aso", "oligonucleotide", "oligonucleotides", "mrna", "crispr",
          "inhibitor", "inhibitors", "agonist", "agonists", "antagonist", "antagonists", "blocker", "blockers",
          "modulator", "modulators", "analogue", "analog", "analogues", "analogs", "mimetic", "mimetics",
          "peptide", "peptides", "supplement", "supplements", "drug", "drugs", "therapy", "therapies",
          "treatment", "treatments", "agent", "agents", "compound", "compounds", "antibody", "antibodies",
          "monoclonal", "statin", "statins", "vaccine", "generation", "new", "novel", "next"}
# doses only: a number with a unit, or a bare number between spaces. Digits glued to letters or after a
# hyphen belong to the compound name (PCSK9, BPC-157, MK-7, SS-31, vitamin K2) and must survive.
_DOSE_RE = re.compile(r"(?<![A-Za-z0-9-])\d+(?:[.,]\d+)?\s*[-–]\s*\d+(?:[.,]\d+)?\s*(?:mg|g|kg|mcg|µg|ug|iu|ui|ml|l|%)?"
                      r"|(?<![A-Za-z0-9-])\d+(?:[.,]\d+)?\s*(?:mg|g|kg|mcg|µg|ug|iu|ui|ml|l|%)?(?![A-Za-z0-9-])", re.I)


def query_terms(name: str, limit: int = 3) -> list[str]:
    """Turn a free-form intervention name into 1-3 searchable compound terms.

    The agent writes things like "Rosuvastatin low-dose", "PCSK9 inhibitors / inclisiran" or
    "Cholecalciferol (vitamin D3) plus magnesium glycinate": a literal TITLE search finds nothing,
    which used to grade well-established drugs as "no evidence".
    """
    name = re.sub(r"\((.*?)\)", r" , \1 , ", name)      # parentheses become alternatives
    out: list[str] = []
    for part in _SPLIT_RE.split(name):
        part = _QUALIFIER_RE.sub(" ", part)
        part = _DOSE_RE.sub(" ", part)
        part = re.sub(r"[^A-Za-z0-9\-\s]", " ", part)
        part = re.sub(r"(?<![A-Za-z0-9])-|-(?![A-Za-z0-9])", " ", part)   # dangling hyphens
        part = " ".join(part.split()).strip(" -")
        if not (3 <= len(part) <= 45):
            continue
        if part.lower() in _GENERIC or part.lower() in (o.lower() for o in out):
            continue
        tokens = [t for t in re.split(r"[\s-]+", part.lower()) if t]
        if all(t in _CLASS or t in _GENERIC for t in tokens):
            continue   # pure modality/class term: not a searchable compound
        out.append(part)
        if len(out) >= limit:
            break
    return out   # empty means "no searchable compound in this name" (a modality, a habit, a protocol)


def _clean(name: str) -> str:
    terms = query_terms(name, limit=1)
    return terms[0] if terms else ""


async def _epmc(client: httpx.Client | httpx.AsyncClient, query: str, page_size: int = 5) -> tuple[int, list[dict]]:
    r = await client.get(EPMC, params={"query": query, "format": "json", "pageSize": page_size, "resultType": "lite", "sort": "CITED desc"}, timeout=40)
    r.raise_for_status()
    d = r.json()
    return int(d.get("hitCount", 0)), d.get("resultList", {}).get("result", [])


async def audit_compound(client: httpx.AsyncClient, name: str, sem: asyncio.Semaphore) -> EvidenceAudit:
    """Query Europe PMC / ClinicalTrials.gov for each candidate term and keep the best-supported one."""
    terms = query_terms(name)
    if not terms:
        return EvidenceAudit(compound=name, error="no evaluable (clase/modalidad o hábito, no un compuesto concreto)")
    best: EvidenceAudit | None = None
    for term in terms:
        a = await _audit_term(client, name, term, sem)
        if best is None or (a.rct_meta_count, a.human_recent_count, a.total_count) > \
                           (best.rct_meta_count, best.human_recent_count, best.total_count):
            best = a
        if a.rct_meta_count >= 20:      # already an A, no need to try the rest
            break
    return best or EvidenceAudit(compound=name, error="sin término de búsqueda")


async def _audit_term(client: httpx.AsyncClient, name: str, term: str, sem: asyncio.Semaphore) -> EvidenceAudit:
    a = EvidenceAudit(compound=name, term=term)
    if len(term) < 3:
        a.error = "nombre demasiado corto"
        return a
    year = date.today().year
    title_q = f'(TITLE:"{term}")'
    async with sem:
        try:
            a.total_count, top = await _epmc(client, f'{title_q} AND SRC:MED', 5)
            a.rct_meta_count, top_rct = await _epmc(
                client, f'{title_q} AND SRC:MED AND (PUB_TYPE:"Randomized Controlled Trial" OR PUB_TYPE:"Meta-Analysis" OR PUB_TYPE:"Systematic Review")', 5)
            a.human_recent_count, _ = await _epmc(
                client, f'{title_q} AND SRC:MED AND (KW:human OR MESH:"Humans") AND FIRST_PDATE:[{year-5}-01-01 TO {year}-12-31]', 1)
            a.top_titles = [f"{t.get('title','').rstrip('.')} ({t.get('journalTitle','')}, {t.get('pubYear','')}; PMID {t.get('pmid','')})"
                            for t in (top_rct or top)[:4]]
        except Exception as e:  # noqa: BLE001
            a.error = f"EuropePMC: {e!s:.80}"
        try:
            r = await client.get(CTG, params={"query.intr": term, "filter.overallStatus": "RECRUITING,ACTIVE_NOT_RECRUITING,ENROLLING_BY_INVITATION,NOT_YET_RECRUITING",
                                              "pageSize": 4, "countTotal": "true", "fields": "NCTId,BriefTitle,Phase,OverallStatus"}, timeout=40)
            if r.status_code == 200:
                d = r.json()
                a.active_trials = int(d.get("totalCount", 0))
                for st in d.get("studies", []):
                    idm = st.get("protocolSection", {}).get("identificationModule", {})
                    dm = st.get("protocolSection", {}).get("designModule", {})
                    a.trial_titles.append(f"{idm.get('nctId')}: {idm.get('briefTitle','')[:90]} [{', '.join(dm.get('phases', []) or [])}]")
        except Exception as e:  # noqa: BLE001
            a.error = (a.error + "; " if a.error else "") + f"CT.gov: {e!s:.60}"
    return a


async def audit_many(names: list[str], concurrency: int = 4) -> dict[str, EvidenceAudit]:
    names = list(dict.fromkeys(n for n in names if n))
    sem = asyncio.Semaphore(concurrency)
    async with httpx.AsyncClient(headers={"User-Agent": "genostack/0.1"}, timeout=40) as client:
        res = await asyncio.gather(*(audit_compound(client, n, sem) for n in names))
    return {a.compound: a for a in res}


async def ensembl_frequencies(rsids: list[str]) -> dict[str, dict]:
    """Global + EUR allele frequencies (1000 Genomes) for a handful of rsIDs via Ensembl REST (GRCh37 server)."""
    out: dict[str, dict] = {}
    rsids = [r for r in dict.fromkeys(rsids) if r.startswith("rs")][:200]
    if not rsids:
        return out
    async with httpx.AsyncClient(timeout=60, headers={"Content-Type": "application/json", "Accept": "application/json"}) as client:
        for i in range(0, len(rsids), 200):
            chunk = rsids[i:i + 200]
            try:
                r = await client.post("https://grch37.rest.ensembl.org/variation/homo_sapiens?pops=1", json={"ids": chunk})
                if r.status_code != 200:
                    continue
                for rid, v in r.json().items():
                    freqs = {}
                    for p in v.get("populations", []):
                        if p.get("population") in ("1000GENOMES:phase_3:ALL", "1000GENOMES:phase_3:EUR", "gnomADe:ALL"):
                            freqs.setdefault(p["population"].split(":")[-1] if "1000" in p["population"] else "gnomAD", {})[p["allele"]] = p["frequency"]
                    out[rid] = {"maf": v.get("MAF"), "minor_allele": v.get("minor_allele"), "ancestral": v.get("ancestral_allele"),
                                "freqs": freqs, "consequence": v.get("most_severe_consequence")}
            except Exception:  # noqa: BLE001
                continue
    return out
