"""Core tests: parsers, strand handling, panel rules, ranking, ACP client with a mock agent."""
from __future__ import annotations

import asyncio
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genostack.parsers import Genome, Call, parse_raw, complement  # noqa: E402
from genostack.panel import load_panel, HaploRule, _eval_cond  # noqa: E402
from genostack.annotate import Finding, eval_panel, eval_haplotypes  # noqa: E402
from genostack.score import rank  # noqa: E402
from genostack.acp_client import AgentRunner, extract_json  # noqa: E402

FIXTURE = ROOT / "tests" / "fixture_23andme.txt"


def test_parse_23andme_fixture():
    g = parse_raw(FIXTURE)
    assert g.source_format == "23andMe"
    assert g.build == "GRCh37"
    assert g.get("rs1801133") == "AG"
    assert g.get("rs6323") == "T"          # haploid X call kept as single letter
    assert g.n_nocalls > 0


def test_parse_ancestry_and_csv(tmp_path):
    a = tmp_path / "anc.txt"
    a.write_text("#AncestryDNA raw data\nrsid\tchromosome\tposition\tallele1\tallele2\nrs1801133\t1\t11856378\tG\tA\nrs2\t23\t100\t0\t0\n", encoding="utf-8")
    g = parse_raw(a)
    assert g.source_format == "AncestryDNA" and g.get("rs1801133") == "AG" and g.get("rs2") is None
    m = tmp_path / "mh.csv"
    m.write_text('##fileformat=MyHeritage\nRSID,CHROMOSOME,POSITION,RESULT\n"rs1801133","1","11856378","GA"\n', encoding="utf-8")
    g2 = parse_raw(m)
    assert g2.get("rs1801133") == "AG"


def test_complement():
    assert complement("AG") == "CT"
    assert complement("A") == "T"


def test_panel_loads_and_is_consistent():
    p = load_panel()
    assert len(p.variants) >= 40
    ids = [v.id for v in p.variants]
    assert len(ids) == len(set(ids))
    for v in p.variants:
        assert v.effect_allele in "ACGTDI" and len(v.effect_allele) == 1, v.id
        assert v.other_allele in "ACGTDI" and len(v.other_allele) == 1, v.id
        # all three diploid genotypes present
        a, b = sorted([v.effect_allele, v.other_allele])
        for gt in (a + a, a + b, b + b):
            assert gt in v.genotypes, f"{v.id} falta genotipo {gt}"
        for gt, g in v.genotypes.items():
            assert 0 <= float(g["impact"]) <= 4, v.id
    for r in p.rules:
        assert r.type in ("combo", "allele_count"), r.id
        assert r.rsids, r.id


def test_panel_eval_on_fixture():
    p = load_panel()
    g = parse_raw(FIXTURE)
    fs = eval_panel(g, p.variants)
    assert any(f.rsid == "rs1801133" for f in fs)  # MTHFR het has impact >= 1
    hs = eval_haplotypes(g, p.rules)
    names = {f.extra["rule_id"]: f.genotype for f in hs}
    # fixture: rs429358 CT + rs7412 CC -> APOE e3/e4
    if "APOE" in names:
        assert "4" in names["APOE"]


def test_allele_count_rule():
    rule = HaploRule(id="X", module="pharmaco", gene="CYP2C19", title="t", type="allele_count", rsids=["rs4244285", "rs12248560"],
                     combos={}, alleles=[{"rsid": "rs4244285", "allele": "A", "function": "loss"}, {"rsid": "rs12248560", "allele": "T", "function": "gain"}],
                     phenotypes=[{"when": "loss>=2", "name": "PM", "impact": 4}, {"when": "loss==1", "name": "IM", "impact": 2},
                                 {"when": "gain>=1", "name": "RM", "impact": 2}, {"when": "true", "name": "NM", "impact": 0}],
                     evidence="A", mechanism="", levers=[])
    g = Genome("t", "GRCh37", {"rs4244285": Call("10", 1, "AG"), "rs12248560": Call("10", 2, "CC")}, 2, 0)
    assert rule.evaluate(g)["name"] == "IM"
    g = Genome("t", "GRCh37", {"rs4244285": Call("10", 1, "AA")}, 1, 0)
    r = rule.evaluate(g)
    assert r["name"] == "PM" and r["partial"] is True
    assert _eval_cond("loss>=1 and gain==0", {"loss": 1, "gain": 0, "count": 1})


def test_rank_groups_modules():
    p = load_panel()
    g = parse_raw(FIXTURE)
    fs = eval_panel(g, p.variants) + eval_haplotypes(g, p.rules)
    fs.append(Finding(source="gwas", rsid="rs999", gene="ZZZ9", genotype="AA", title="Some trait", detail="", score=6.0, evidence="B"))
    an = rank(fs, p, ["health"])
    assert any(b.id == "other" for b in an.modules)
    assert all(f.extra.get("final") is not None for b in an.modules for f in b.all_findings())


def test_extract_json():
    assert extract_json("bla ```json\n{\"a\": 1,}\n``` fin")["a"] == 1
    assert extract_json("{\"b\": [1,2]}")["b"] == [1, 2]
    with pytest.raises(ValueError):
        extract_json("nada")


def test_acp_mock_agent_roundtrip():
    async def main():
        cmd = f'"{sys.executable}" "{ROOT / "tests" / "mock_agent.py"}"'
        async with AgentRunner(cmd, progress=lambda l, t: None, concurrency=2) as ar:
            assert "mock" in ar.agent_info
            texts = await asyncio.gather(*(ar.ask('{"module": "m%d"}' % i, label=f"m{i}") for i in range(3)))
            for i, t in enumerate(texts):
                assert extract_json(t)["module"] == f"m{i}"
    asyncio.run(main())


def test_annotate_with_local_db():
    from genostack.config import DB_PATH
    if not DB_PATH.exists():
        pytest.skip("base local no construida (genostack fetch-data)")
    from genostack.annotate import annotate, open_db
    g = parse_raw(FIXTURE)
    g.calls["i3002432"] = Call("11", 46761055, "AG")  # 23andMe internal id for rs1799963 (F2 G20210A)
    con = open_db()
    try:
        findings, stats = annotate(g, con)
    finally:
        con.close()
    assert stats["n_in_gwas"] > 0 and stats["n_in_clinvar"] > 0
    assert any(f.source == "haplotype" for f in findings)
    if stats.get("snpedia_available"):
        assert any(f.source == "snpedia" for f in findings)
        assert g.get("rs1799963") == "AG"  # alias resolved through snpedia_alias


def test_query_terms_extracts_real_compounds():
    """The agent writes phrases; the audit must search the actual compound (and refuse pure class words)."""
    from genostack.evidence import query_terms
    assert query_terms("Rosuvastatin low-dose") == ["Rosuvastatin"]
    assert query_terms("PCSK9 inhibitors / inclisiran") == ["PCSK9 inhibitors", "inclisiran"]
    assert query_terms("Creatine monohydrate 5 g/day") == ["Creatine"]
    assert query_terms("Urolithin A 500-1000 mg") == ["Urolithin A"]
    # digits that belong to the compound name must survive
    for name in ("BPC-157", "Menaquinone-7", "SS-31"):
        assert name in query_terms(name)[0]
    # a modality is not a compound: searching it would grade an experimental drug as established
    assert query_terms("siRNA") == []
    assert query_terms("Novel peptide therapy") == []


def test_transient_and_quota_classification():
    from genostack.acp_client import QUOTA_RE, TRANSIENT_RE
    for msg in ("API Error: 529 Overloaded", "The response stopped arriving",
                "El agente terminó con stop_reason=end_turn sin texto (sesión x)"):
        assert TRANSIENT_RE.search(msg), msg
    assert QUOTA_RE.search("usage limit reached")
    assert not TRANSIENT_RE.search("400 invalid request")
