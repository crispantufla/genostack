"""Render the final report as Markdown (printed to console via rich, and saved next to the input)."""
from __future__ import annotations

from datetime import datetime

from rich.console import Console
from rich.markdown import Markdown

from .annotate import Finding
from .research import ModuleResearch, Synthesis
from .score import Analysis

DIR_ICON = {"risk": "⚠", "beneficial": "✅", "mixed": "◐", "info": "ℹ"}
SRC_LABEL = {"panel": "panel", "haplotype": "haplotipo", "snpedia": "SNPedia", "clinvar": "ClinVar", "gwas": "GWAS", "pgx": "PharmGKB"}


def _esc(s: str) -> str:
    return (s or "").replace("|", "/").replace("\n", " ").strip()


def _finding_line(f: Finding) -> str:
    icon = DIR_ICON.get(f.direction, "•")
    if f.source == "haplotype":
        gts = f.extra.get("genotypes") or {}
        head = f"{icon} **{_esc(f.title)}** — " + " ".join(f"`{r}`={g}" for r, g in gts.items())
    else:
        head = f"{icon} **{_esc(f.title)}** — `{f.rsid}` {f.genotype}"
    if f.partial:
        head += " _(haplotipo parcial)_"
    body = f"  \n  {_esc(f.detail)}"
    meta = f"  \n  _impacto {f.impact:.1f}/4 · evidencia {f.evidence} · {SRC_LABEL.get(f.source, f.source)}_"
    if f.effect_size:
        meta += f" · efecto: {_esc(f.effect_size)}"
    out = head + body + meta
    if f.extra.get("freq") and f.extra["freq"].get("maf") is not None:
        out += f" · MAF global {f.extra['freq']['maf']}"
    if f.levers:
        out += "  \n  Palancas semilla: " + "; ".join(_esc(x) for x in f.levers[:6])
    if f.caveats and f.source in ("clinvar", "haplotype"):
        out += f"  \n  _{_esc(f.caveats)}_"
    return out


def _table(rows: list[list[str]], header: list[str]) -> str:
    out = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for r in rows:
        out.append("| " + " | ".join(_esc(str(c)) for c in r) + " |")
    return "\n".join(out)


def render_markdown(analysis: Analysis, research: list[ModuleResearch] | None, synth: Synthesis | None, stats: dict,
                    meta: dict, lang: str = "es") -> str:
    md: list[str] = []
    md.append(f"# GENOSTACK — Informe genético y de intervenciones\n")
    md.append(f"_Generado {datetime.now():%Y-%m-%d %H:%M} · archivo: {meta.get('file')} · formato {meta.get('format')} · build {meta.get('build')} · sexo inferido: {meta.get('sex')}_\n")
    md.append(f"**Cobertura:** {stats.get('n_calls', 0):,} SNPs leídos · {stats.get('n_in_snpedia', 0):,} en SNPedia · {stats.get('n_in_clinvar', 0):,} en ClinVar · "
              f"{stats.get('n_in_gwas', 0):,} en GWAS Catalog · {stats.get('n_in_pgx', 0):,} en PharmGKB · panel curado {stats.get('panel_covered', 0)}/{stats.get('panel_variants', 0)} variantes + {stats.get('panel_rules', 0)} reglas de haplotipo.  \n"
              f"**Objetivos ponderados:** {', '.join(analysis.goals)} · **Agente:** {meta.get('agent') or 'ninguno (modo offline)'}\n")
    if not stats.get("snpedia_available"):
        md.append("> ⚠ SNPedia no está en la base local (ejecuta `genostack fetch-data` completo para añadirla).\n")
    md.append("> ⚠ **Limitaciones:** arrays de consumo (sin imputación, alta tasa de falsos positivos en variantes raras → confirmar ClinVar patogénicas con secuenciación clínica; sin CNV como GSTM1/GSTT1 null; sin tipado HLA real; SNPs A/T-C/G con orientación incierta se marcan). Esto no es consejo médico: discútelo con un profesional antes de actuar, especialmente fármacos y compuestos experimentales.\n")

    # ---- executive summary
    if synth and not synth.error:
        md.append("## 1. Resumen ejecutivo\n")
        md.append(synth.executive_summary.strip() + "\n")
        if synth.genetic_profile_summary:
            md.append("**Perfil genético integrado:** " + synth.genetic_profile_summary.strip() + "\n")
        if synth.core_stack:
            md.append("### Stack central\n")
            rows = [[f"{i.priority}", f"**{i.name}**", i.category, i.dose, i.timing, i.evidence_grade, ", ".join(i.genes[:4]), i.why] for i in sorted(synth.core_stack, key=lambda x: x.priority)]
            md.append(_table(rows, ["Prio", "Intervención", "Tipo", "Dosis", "Momento", "Evid.", "Genes", "Por qué (vínculo genético)"]) + "\n")
        if synth.advanced_experimental:
            md.append("### Avanzado / experimental (sólo con supervisión; revisar estatus legal)\n")
            rows = [[f"{i.priority}", f"**{i.name}**", i.category, i.dose, i.evidence_grade, ", ".join(i.genes[:4]), i.why] for i in sorted(synth.advanced_experimental, key=lambda x: x.priority)]
            md.append(_table(rows, ["Prio", "Intervención", "Tipo", "Dosis", "Evid.", "Genes", "Por qué / estatus"]) + "\n")
    elif synth and synth.error:
        md.append(f"## 1. Resumen ejecutivo\n\n> La síntesis del agente falló: {synth.error}. Se muestran los módulos individuales.\n")
    else:
        md.append("## 1. Resumen ejecutivo (modo offline)\n")
        top = sorted((f for b in analysis.modules for f in b.findings),
                     key=lambda f: (f.source not in ("panel", "haplotype", "clinvar"), -f.extra.get("final", 0)))[:15]
        for f in top:
            gt = "" if f.source == "haplotype" else f" ({f.genotype})"
            md.append(f"- {DIR_ICON.get(f.direction,'•')} **{_esc(f.title)}**{gt} — {_esc(f.detail)[:160]}")
        md.append("")

    # ---- modules
    md.append("## 2. Módulos funcionales\n")
    rmap = {r.module: r for r in (research or [])}
    for i, b in enumerate(analysis.modules, 1):
        if b.id == "other":
            continue
        md.append(f"### 2.{i} {b.title}  \n_peso {b.weight} · genes: {', '.join(b.genes[:14])}_\n")
        if b.description:
            md.append(f"_{_esc(b.description)}_\n")
        r = rmap.get(b.id)
        if r and r.interpretation:
            md.append("**Interpretación:** " + r.interpretation.strip() + "\n")
        md.append("**Hallazgos principales:**\n")
        for f in b.findings:
            md.append("- " + _finding_line(f))
        if b.minor:
            md.append("\n<details><summary>Variantes de efecto leve (contexto)</summary>\n")
            for f in b.minor:
                md.append(f"- {_esc(f.title)} `{f.rsid}` {f.genotype}: {_esc(f.detail)[:140]}")
            md.append("\n</details>\n")
        if r and r.interventions:
            md.append("\n**Intervenciones (investigación del agente, ordenadas por prioridad):**\n")
            rows = []
            for it in sorted(r.interventions, key=lambda x: x.priority):
                aud = it.audit or {}
                audit_txt = (f"«{aud.get('term','?')}»: RCT/meta {aud.get('rct_meta','?')} · hum.5a {aud.get('human_recent','?')} · "
                             f"ensayos {aud.get('active_trials','?')} → {aud.get('grade_hint','?')}") if aud else "—"
                rows.append([str(it.priority), f"**{it.name}**", it.category, ", ".join(it.target_genes[:3]), it.dose, it.effect_magnitude or it.expected_effect,
                             f"{it.evidence_grade} ({it.novelty})", audit_txt, it.risks_interactions[:160], "; ".join(it.key_studies[:2])])
            md.append(_table(rows, ["Prio", "Intervención", "Tipo", "Genes", "Dosis", "Efecto esperado", "Evid.", "Auditoría EPMC/CT.gov", "Riesgos/interacciones", "Estudios"]) + "\n")
            for it in sorted(r.interventions, key=lambda x: x.priority)[:6]:
                md.append(f"- **{it.name}** — {_esc(it.mechanism)[:300]}{' · Estatus: ' + _esc(it.legal_status) if it.legal_status else ''}")
            md.append("")
            if r.avoid:
                md.append("**Evitar / cuidado por genotipo:** " + "; ".join(_esc(x) for x in r.avoid) + "\n")
            if r.labs:
                md.append("**Analíticas:** " + "; ".join(_esc(x) for x in r.labs) + "\n")
            if r.open_questions:
                md.append("_Preguntas abiertas:_ " + "; ".join(_esc(x) for x in r.open_questions) + "\n")
        elif r and r.error:
            md.append(f"> Investigación del agente no disponible para este módulo: {r.error}\n")
        md.append("")

    # ---- other findings
    other = [b for b in analysis.modules if b.id == "other"]
    if other and other[0].findings:
        md.append("## 3. Otros hallazgos significativos (cruce automático fuera del panel)\n")
        md.append("_SNPedia magnitud ≥ 2.5, ClinVar ≥ 2★, GWAS p ≤ 5e-8 con OR ≥ 1.3 o |β| ≥ 0.1 SD, PharmGKB 1A-2B. Revisar con criterio: son asociaciones, no diagnósticos._\n")
        for f in other[0].findings:
            md.append("- " + _finding_line(f))
        md.append("")

    # ---- pharmaco table
    ph = [b for b in analysis.modules if b.id == "pharmaco"]
    if ph:
        md.append("## 4. Farmacogenómica — tabla rápida\n")
        rows = [[f.gene, f.genotype, _esc(f.title), _esc(f.detail)[:140], f.evidence] for f in ph[0].all_findings()]
        md.append(_table(rows, ["Gen", "Genotipo/diplotipo", "Fenotipo", "Detalle", "Evid."]) + "\n")

    # ---- synthesis details
    if synth and not synth.error:
        md.append("## 5. Integración: conflictos, orden, analíticas\n")
        if synth.conflicts:
            md.append("**Conflictos e interacciones:**\n" + "\n".join(f"- {_esc(x)}" for x in synth.conflicts) + "\n")
        if synth.avoid:
            md.append("**Evitar / limitar:**\n" + "\n".join(f"- {_esc(x)}" for x in synth.avoid) + "\n")
        if synth.introduction_order:
            md.append("**Orden de introducción:**\n" + "\n".join(f"{i+1}. {_esc(x)}" for i, x in enumerate(synth.introduction_order)) + "\n")
        if synth.labs:
            md.append("**Analíticas (basal + seguimiento):**\n" + "\n".join(f"- {_esc(x)}" for x in synth.labs) + "\n")
        if synth.confirm_with_sequencing:
            md.append("**Confirmar con secuenciación clínica:**\n" + "\n".join(f"- {_esc(x)}" for x in synth.confirm_with_sequencing) + "\n")
        if synth.unexpected_findings:
            md.append("**Hallazgos inesperados a vigilar:**\n" + "\n".join(f"- {_esc(x)}" for x in synth.unexpected_findings) + "\n")
        if synth.caveats:
            md.append("_Limitaciones:_ " + _esc(synth.caveats) + "\n")

    md.append("## 6. Fuentes\n")
    md.append(f"Bases locales: {meta.get('db_built', '?')} — SNPedia (CC-BY-NC-SA), ClinVar, GWAS Catalog, Open Targets Platform (incl. PharmGKB/ClinPGx PGx), DGIdb, CPIC. "
              "Online: Europe PMC, ClinicalTrials.gov, Ensembl REST y las búsquedas del agente ACP. "
              f"Artefactos del agente (prompts/respuestas): {meta.get('agent_log_dir') or '—'}\n")
    return "\n".join(md)


def print_report(md: str, console: Console | None = None) -> None:
    console = console or Console()
    console.print(Markdown(md))
