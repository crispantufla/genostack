"""Research orchestration: per-module prompts to the ACP agent, deterministic evidence audit, final synthesis."""
from __future__ import annotations

import asyncio
import json
from typing import Any

from pydantic import BaseModel, Field, ValidationError

from .acp_client import AgentAuthError, AgentQuotaError, AgentRunner, extract_json
from .annotate import Finding
from .evidence import EvidenceAudit, audit_many
from .score import Analysis, ModuleBlock


# ----------------------------------------------------------------------------- schemas
class Intervention(BaseModel):
    name: str
    name_en: str = ""                     # English/INN name used for the literature audit
    category: str = "suplemento"          # suplemento | péptido | fármaco_rx | experimental | estilo_de_vida | dieta
    target_genes: list[str] = Field(default_factory=list)
    rsids: list[str] = Field(default_factory=list)
    mechanism: str = ""
    expected_effect: str = ""
    effect_magnitude: str = ""            # cuantificado cuando sea posible
    evidence_grade: str = "C"             # A-E
    key_studies: list[str] = Field(default_factory=list)
    dose: str = ""
    risks_interactions: str = ""
    novelty: str = "establecido"          # establecido | novedoso | experimental
    legal_status: str = ""
    priority: int = 3                     # 1 (máxima) - 5
    audit: dict | None = None             # filled by evidence audit


class ModuleResearch(BaseModel):
    module: str
    interpretation: str = ""
    interventions: list[Intervention] = Field(default_factory=list)
    avoid: list[str] = Field(default_factory=list)
    labs: list[str] = Field(default_factory=list)
    open_questions: list[str] = Field(default_factory=list)
    error: str | None = None


class StackItem(BaseModel):
    name: str
    why: str = ""
    dose: str = ""
    timing: str = ""
    evidence_grade: str = "C"
    priority: int = 3
    category: str = ""
    genes: list[str] = Field(default_factory=list)


class Synthesis(BaseModel):
    executive_summary: str = ""
    genetic_profile_summary: str = ""
    core_stack: list[StackItem] = Field(default_factory=list)
    advanced_experimental: list[StackItem] = Field(default_factory=list)
    avoid: list[str] = Field(default_factory=list)
    conflicts: list[str] = Field(default_factory=list)
    introduction_order: list[str] = Field(default_factory=list)
    labs: list[str] = Field(default_factory=list)
    confirm_with_sequencing: list[str] = Field(default_factory=list)
    unexpected_findings: list[str] = Field(default_factory=list)
    caveats: str = ""
    error: str | None = None


# ----------------------------------------------------------------------------- prompt building
LANG = {
    "es": {"lang_name": "español"},
    "en": {"lang_name": "English"},
}


def _finding_row(f: Finding) -> dict:
    row = {
        "rsid": f.rsid, "gene": f.gene, "genotype": f.genotype, "title": f.title, "phenotype": f.detail[:350],
        "impact_0_4": round(f.impact, 1), "evidence": f.evidence, "direction": f.direction, "source": f.source,
    }
    if f.effect_size:
        row["effect_size"] = f.effect_size
    if f.mechanism:
        row["mechanism"] = f.mechanism[:300]
    if f.levers:
        row["seed_levers"] = f.levers[:8]
    if f.caveats:
        row["caveats"] = f.caveats[:200]
    if f.partial:
        row["partial_haplotype"] = True
    if f.extra.get("support"):
        row["db_support"] = f.extra["support"][:5]
    if f.extra.get("freq"):
        row["population_freq"] = f.extra["freq"]
    return row


def _gene_ctx_compact(ctx: dict[str, dict]) -> dict:
    out = {}
    for g, c in ctx.items():
        if not c:
            continue
        item = {}
        if c.get("function"):
            item["function"] = c["function"][:500]
        if c.get("ot_drugs"):
            item["known_or_trial_drugs_for_target"] = [f"{d['drug']} [{d.get('max_stage') or '?'}; {d.get('action') or ''} {d.get('mechanism') or ''}]".strip() for d in c["ot_drugs"][:15]]
        if c.get("ot_moa"):
            extra = [f"{d['drug']} [{d.get('max_stage') or '?'}; {d.get('action') or ''}]" for d in c["ot_moa"][:15]
                     if not any(d["drug"] == x["drug"] for x in c.get("ot_drugs", []))]
            if extra:
                item["other_molecules_with_moa_on_target"] = extra
        if c.get("dgidb"):
            item["dgidb_interactions"] = [f"{d['drug']}{' ('+d['type']+')' if d.get('type') else ''}{' ✓aprobado' if d.get('approved') else ''}" for d in c["dgidb"][:15]]
        if c.get("cpic"):
            item["cpic_guidelines"] = [f"{d['drug']} (nivel {d['level']})" for d in c["cpic"][:12]]
        if c.get("pathways"):
            item["pathways"] = c["pathways"][:6]
        if item:
            out[g] = item
    return out


def module_prompt(block: ModuleBlock, gene_ctx: dict[str, dict], profile_summary: str, lang: str, goals: list[str]) -> str:
    L = LANG.get(lang, LANG["es"])
    payload = {
        "module": block.id, "title": block.title, "module_goals": block.goals, "user_goals": goals,
        "findings_major": [_finding_row(f) for f in block.findings],
        "findings_minor_context": [_finding_row(f) for f in block.minor],
        "gene_context_from_local_db": _gene_ctx_compact(gene_ctx),
        "whole_genome_profile_summary": profile_summary,
    }
    schema = {
        "module": block.id,
        "interpretation": "4-8 frases: qué dice este módulo del perfil de ESTA persona (fenotipo integrado, no lista de SNPs)",
        "interventions": [{
            "name": "nombre del compuesto/intervención", "name_en": "SOLO el nombre del principio activo en inglés/INN, 1-3 palabras, sin dosis ni adjetivos (p.ej. 'rosuvastatin', 'inclisiran', 'menaquinone-7'; NO 'low-dose rosuvastatin as first-line statin')", "category": "suplemento | péptido | fármaco_rx | experimental | estilo_de_vida | dieta",
            "target_genes": ["GEN"], "rsids": ["rs..."], "mechanism": "cómo corrige/compensa/amplifica la variante",
            "expected_effect": "efecto esperado concreto", "effect_magnitude": "cuantifica: p.ej. 'homocisteína −25 %', 'SMD 0.4', 'OR 0.7'",
            "evidence_grade": "A|B|C|D|E", "key_studies": ["Autor año, revista (PMID/DOI/NCT)"], "dose": "dosis, forma, momento",
            "risks_interactions": "riesgos, interacciones, contraindicaciones", "novelty": "establecido | novedoso | experimental",
            "legal_status": "OTC / receta / investigacional / gris", "priority": 1,
        }],
        "avoid": ["qué evitar o dosificar con cuidado POR ESTE genotipo"],
        "labs": ["analíticas/biomarcadores que confirman o guían"],
        "open_questions": ["incertidumbres relevantes"],
    }
    return f"""Eres un experto en nutrigenómica, farmacogenómica, farmacología de péptidos, medicina de longevidad y optimización del rendimiento (energía, cognición, productividad). Investigas intervenciones para UNA persona concreta a partir de su genotipo. Responde en {L['lang_name']}.

## Perfil genético de este módulo ({block.title})
```json
{json.dumps(payload, ensure_ascii=False, indent=1)}
```
Los `seed_levers` son sólo un punto de partida curado: verifícalos, matízalos y amplíalos. `gene_context_from_local_db` lista fármacos/moléculas conocidas y en ensayo para las dianas (Open Targets/DGIdb/CPIC) como inspiración para "amplificar o modular" el gen; no todas son relevantes ni seguras.

## Tarea
Haz una investigación EXHAUSTIVA (usa tus herramientas de búsqueda web y lectura de páginas tanto como haga falta: PubMed/Europe PMC, revisiones 2023-2026, preprints bioRxiv/medRxiv, ClinicalTrials.gov, Examine, guías CPIC/PharmGKB) y devuelve las intervenciones con EFECTO SIGNIFICATIVO para ESTE perfil:
1. Cubre obligatoriamente y por separado: (a) suplementos/nutracéuticos, (b) péptidos (incl. los de uso en optimización/longevidad si hay base mecanística o datos humanos: p.ej. análogos GLP-1, MOTS-c, SS-31/elamipretida, humanina, BPC-157, epitalón, semax/selank, GHK-Cu, tesamorelina, etc. SÓLO si encajan con el genotipo), (c) fármacos con receta (incl. uso off-label razonado: metformina, rapamicina, acarbosa, SGLT2i, litio microdosis, estatinas/PCSK9i, etc. SÓLO si encajan), (d) compuestos experimentales/novedosos (2023-2026: nuevas moléculas, ensayos en curso, formulaciones nuevas), (e) dieta/estilo de vida de alto impacto.
2. Filtro de significancia: incluye sólo lo que tiene efecto clínicamente relevante (cambio de biomarcador u outcome relevante; no "estadísticamente significativo pero marginal"). Para lo experimental exige mecanismo directo sobre el gen/ruta afectada + al menos datos humanos preliminares o animales robustos, y etiquétalo como tal.
3. Cada intervención debe vincularse EXPLÍCITAMENTE a genotipos concretos de arriba (rsids). Nada de "es bueno en general".
4. Sé concreto en dosis, forma, momento, duración y en los riesgos/interacciones de ESTE genotipo (p.ej. si hay COMT lento, avisa de donantes de metilo; si HFE, de hierro/vit C; etc.).
5. Da 6-14 intervenciones, ordenadas por prioridad (1 = máxima). Prioriza por (magnitud de efecto × evidencia × relevancia a los objetivos del usuario: {', '.join(goals)}).
6. Para cada una, cita 1-3 estudios clave con PMID/DOI/NCT reales (si no recuerdas el identificador exacto, da autor+año+revista; no inventes identificadores).

## Formato de salida
Devuelve ÚNICAMENTE un bloque ```json con esta estructura (sin texto antes ni después):
```json
{json.dumps(schema, ensure_ascii=False, indent=1)}
```"""


def synthesis_prompt(analysis: Analysis, modules: list[ModuleResearch], audits: dict[str, EvidenceAudit], profile_summary: str,
                     lang: str, goals: list[str]) -> str:
    L = LANG.get(lang, LANG["es"])
    compact = []
    for m in modules:
        if m.error and not m.interventions:
            compact.append({"module": m.module, "error": m.error})
            continue
        compact.append({
            "module": m.module, "interpretation": m.interpretation,
            "interventions": [{
                "name": i.name, "category": i.category, "genes": i.target_genes, "rsids": i.rsids[:6], "effect": i.expected_effect,
                "magnitude": i.effect_magnitude, "grade": i.evidence_grade, "dose": i.dose, "risks": i.risks_interactions[:200],
                "novelty": i.novelty, "priority": i.priority,
                "evidence_audit": i.audit,
            } for i in m.interventions],
            "avoid": m.avoid, "labs": m.labs,
        })
    other = [b for b in analysis.modules if b.id == "other"]
    other_rows = [_finding_row(f) for b in other for f in b.findings[:20]]
    schema = {
        "executive_summary": "10-15 frases en prosa: lo que más importa de este genoma y el stack central",
        "genetic_profile_summary": "síntesis integrada del perfil (fenotipos clave por sistema)",
        "core_stack": [{"name": "", "why": "vinculado a genes/rsids", "dose": "", "timing": "", "evidence_grade": "A|B|C|D|E", "priority": 1, "category": "", "genes": [""]}],
        "advanced_experimental": [{"name": "", "why": "", "dose": "", "timing": "", "evidence_grade": "", "priority": 3, "category": "péptido|fármaco_rx|experimental", "genes": [""]}],
        "avoid": ["qué evitar o limitar por genotipo"],
        "conflicts": ["conflictos/interacciones entre módulos o entre intervenciones y cómo resolverlos"],
        "introduction_order": ["fase 1 (semanas 0-4): ...", "fase 2: ...", "fase 3: ..."],
        "labs": ["analíticas basales y de seguimiento, con valores objetivo"],
        "confirm_with_sequencing": ["variantes raras/patogénicas a confirmar clínicamente"],
        "unexpected_findings": ["hallazgos fuera del panel que merecen atención"],
        "caveats": "limitaciones importantes",
    }
    return f"""Eres el médico-investigador que integra los resultados de varios módulos de análisis genético de UNA persona en un plan único, priorizado y coherente. Responde en {L['lang_name']}.

## Perfil genético global
{profile_summary}

## Resultados por módulo (ya investigados; `evidence_audit` = recuento objetivo en Europe PMC/ClinicalTrials.gov: rct_meta = RCT+meta-análisis, human_recent = estudios humanos últimos 5 años, active_trials = ensayos activos, grade_hint = grado sugerido por la auditoría)
```json
{json.dumps(compact, ensure_ascii=False)}
```

## Otros hallazgos significativos fuera del panel (cruce automático SNPedia/ClinVar/GWAS)
```json
{json.dumps(other_rows, ensure_ascii=False)}
```

## Tarea
1. Construye el STACK CENTRAL (8-15 intervenciones máximo) con lo de mayor impacto × evidencia para los objetivos: {', '.join(goals)}. Usa la auditoría de evidencia para rebajar lo que el módulo sobrevaloró (p.ej. rct_meta=0 y human_recent=0 → no puede ser grado A/B) y para destacar lo novedoso con ensayos activos.
2. Lista aparte lo AVANZADO/EXPERIMENTAL (péptidos, fármacos off-label, moléculas en ensayo) con su estatus legal y riesgos, sólo si encaja con el genotipo.
3. Detecta CONFLICTOS e interacciones entre módulos (p.ej. COMT lento × donantes de metilo; HFE × vitamina C/hierro; APOE4 × grasas saturadas/ciertas formas de DHA; CYP2C19 PM × omeprazol/clopidogrel/citalopram; G6PD × vitamina C IV/azul de metileno; SLCO1B1 × estatinas; CYP1A2 lento × cafeína+estimulantes; MTHFR × ácido fólico sintético; ADORA2A × cafeína+ansiedad) y di cómo resolverlos. Consolida duplicados entre módulos.
4. Propón el ORDEN DE INTRODUCCIÓN por fases (para poder atribuir efectos), y las ANALÍTICAS basales y de seguimiento con valores objetivo.
5. Señala qué variantes raras/patogénicas deben CONFIRMARSE con secuenciación clínica antes de actuar, y qué hallazgos inesperados fuera del panel merecen atención (o descártalos razonadamente).
6. Puedes usar búsqueda web puntualmente para resolver dudas concretas (no repitas la investigación).

Devuelve ÚNICAMENTE un bloque ```json con esta estructura:
```json
{json.dumps(schema, ensure_ascii=False, indent=1)}
```"""


# ----------------------------------------------------------------------------- orchestration
async def research_module(runner: AgentRunner, block: ModuleBlock, gene_ctx: dict, profile_summary: str, lang: str,
                          goals: list[str]) -> ModuleResearch:
    prompt = module_prompt(block, gene_ctx, profile_summary, lang, goals)
    last_err = ""
    for attempt in range(2):
        try:
            text = await runner.ask(prompt if attempt == 0 else prompt + "\n\nIMPORTANTE: tu respuesta anterior no contenía JSON válido. Devuelve SOLO el bloque ```json conforme al esquema.",
                                    label=f"{block.id}{'' if attempt == 0 else '_retry'}")
            data = extract_json(text)
            data["module"] = block.id
            return ModuleResearch.model_validate(data)
        except (AgentAuthError, AgentQuotaError):
            raise
        except (ValueError, ValidationError, RuntimeError) as e:
            last_err = f"{type(e).__name__}: {e!s:.200}"
    return ModuleResearch(module=block.id, error=last_err)


async def run_research(runner: AgentRunner, analysis: Analysis, gene_ctx_by_module: dict[str, dict], profile_summary: str,
                       lang: str, goals: list[str], log=print) -> tuple[list[ModuleResearch], dict[str, EvidenceAudit], Synthesis]:
    blocks = [b for b in analysis.modules if b.id != "other"]
    results = await asyncio.gather(*(research_module(runner, b, gene_ctx_by_module.get(b.id, {}), profile_summary, lang, goals) for b in blocks))
    # evidence audit (deterministic)
    names = [(i.name_en or i.name) for m in results for i in m.interventions if i.category not in ("estilo_de_vida", "dieta", "analítica")]
    log(f"Auditoría de evidencia: {len(set(names))} compuestos en Europe PMC / ClinicalTrials.gov…")
    audits = await audit_many(names)
    for m in results:
        for i in m.interventions:
            a = audits.get(i.name_en or i.name)
            if a:
                i.audit = {"term": a.term, "rct_meta": a.rct_meta_count, "human_recent": a.human_recent_count,
                           "total": a.total_count, "active_trials": a.active_trials, "grade_hint": a.grade_hint(),
                           "top": a.top_titles[:2]}
    # synthesis
    prompt = synthesis_prompt(analysis, results, audits, profile_summary, lang, goals)
    synth = Synthesis(error="sin síntesis")
    for attempt in range(2):
        try:
            text = await runner.ask(prompt if attempt == 0 else prompt + "\n\nIMPORTANTE: devuelve SOLO el bloque ```json conforme al esquema.",
                                    label="synthesis" + ("" if attempt == 0 else "_retry"), timeout=1800)
            synth = Synthesis.model_validate(extract_json(text))
            break
        except (AgentAuthError, AgentQuotaError):
            raise
        except (ValueError, ValidationError, RuntimeError) as e:
            synth = Synthesis(error=f"{type(e).__name__}: {e!s:.200}")
    return results, audits, synth
