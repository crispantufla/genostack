# GENOSTACK — Plan y brainstorming

Programa local de consola: recibe un raw report de ADN (23andMe / AncestryDNA / MyHeritage / FTDNA / VCF simple), lo analiza contra bases de datos descargadas de antemano (≤ 2 GB), y delega la investigación abierta (suplementos, péptidos, fármacos, compuestos experimentales/novedosos) a un agente por **ACP (Agent Client Protocol)**. Salida: un informe profundo en consola (y copia `.md` junto al archivo de entrada). Nada más.

Uso final:

```bash
genostack fetch-data                 # una vez: descarga e indexa las bases locales (≤ 2 GB)
genostack analyze genome_23andme.txt # imprime el informe
```

---

## 0. Principios de diseño

1. **Simple.** Un CLI, dos subcomandos, un único fichero DuckDB local, un agente ACP. Sin UI, sin servidor, sin "features".
2. **Significancia, no mera existencia.** Cada hallazgo genético y cada intervención pasa un filtro de umbral explícito (§4). Si no lo supera, no aparece. Mejor 25 cosas que importan que 400 ruidosas.
3. **Genotipo-específico.** Las intervenciones se buscan *para este genoma*, no listas genéricas de longevidad/nootrópicos. El agente recibe el contexto genético y la evidencia local ya filtrada, y se le exige vincular cada recomendación a variantes concretas.
4. **Exhaustivo en lo que importa.** Exhaustividad = (a) cruce sistemático de TODOS los SNPs del archivo contra todas las bases locales, (b) panel curado de variantes accionables de alto impacto, (c) investigación online abierta y profunda de intervenciones, incl. compuestos nuevos (últimos 2-3 años, preprints, ensayos en curso).
5. **Privacidad.** El genoma crudo nunca sale de la máquina. Al agente sólo viajan los rsIDs+genotipos ya priorizados (decenas, no cientos de miles).
6. **Determinista donde se pueda, LLM donde haga falta.** Parseo, cruces, haplotipos, scoring y conteo de evidencia: código. Síntesis, razonamiento mecanístico y descubrimiento de compuestos novedosos: agente.

---

## 1. Stack técnico

- **Python 3.12+** (uv), `duckdb` (almacén local único, joins rápidos sobre ~600k-1M rsIDs), `httpx` (descargas/APIs), `rich` (consola), `pydantic` (esquemas JSON del agente), `pyyaml` (panel curado).
- **ACP**: paquete `agent-client-protocol` (PyPI; `spawn_agent_process`, `initialize`, `new_session`, `prompt`, handlers `session_update` / `request_permission`). Cliente = nuestro programa; agente = subproceso por stdio. Por defecto `npx -y @agentclientprotocol/claude-agent-acp` (verificado: v0.70.0 resuelve en esta máquina). Configurable a Gemini CLI (`gemini --experimental-acp`), `codex-acp`, etc. con `--agent "<cmd>"`.
- Sin dependencias bio pesadas: el raw report ya trae rsid/cromosoma/posición/genotipo; no hay alineamiento ni imputación.

Estructura:

```
genostack/
  cli.py            fetch-data | analyze
  parsers.py        23andMe v3/v4/v5, AncestryDNA, MyHeritage, FTDNA, VCF mínimo → (rsid, chr, pos, gt)
  fetch.py          descarga + construye data/genostack.duckdb (resumible, presupuesto 2 GB)
  annotate.py       joins rsid ↔ SNPedia/ClinVar/GWAS/PharmGKB; normaliza hebra/orientación
  panel/*.yaml      panel curado de variantes y haplotipos accionables (§3.3)
  haplotypes.py     APOE, CYP2C19/2D6/2C9 (alelos estrella por tag-SNPs), MTHFR compuesto, etc.
  score.py          significancia, filtro, agrupación en módulos funcionales
  evidence.py       Europe PMC / ClinicalTrials.gov / Ensembl REST (online, puntual)
  acp_client.py     spawn agente, permisos, streaming, sesiones concurrentes, timeouts, parse JSON
  research.py       prompts por módulo + auditoría de evidencia + síntesis final
  report.py         render consola (rich) + guarda <input>.genostack.md
tests/              fixtures sintéticos, parser, hebra, haplotipos, scoring, agente mock (echo agent del SDK)
```

---

## 2. Datos offline (fase 1: descargar por adelantado, ≤ 2 GB)

| Fuente | Para qué | Tamaño aprox. | Notas |
|---|---|---|---|
| **SNPedia** (API MediaWiki `bots.snpedia.com/api.php`, caché JSON) | Resumen por SNP y por genotipo, **magnitude** (0-10), **repute**, orientación de hebra, gen | ~250-400 MB | CC-BY-NC-SA (uso personal OK). Es la fuente con el concepto de "magnitud" = exactamente el filtro de significancia que pides. Scrape resumible con rate-limit. |
| **ClinVar** `variant_summary.txt.gz` | Variantes patogénicas/probablemente patogénicas, nivel de revisión (estrellas) | ~350 MB gz | Filtrar a filas con rsID. |
| **GWAS Catalog** associations (TSV completo) | rsid → rasgo, p-valor, OR/beta, alelo de riesgo | ~300 MB | Filtrar p ≤ 5e-8; priorizar rasgos de categorías objetivo (niveles de nutrientes, metabolismo, cognición, sueño, inflamación, cardiometabólico, longevidad). |
| **PharmGKB** (clinicalAnnotations, variantAnnotations, drugLabels, clinicalVariants) | Farmacogenómica: variante → fármaco, nivel 1A-4 | ~50 MB | Base para "fármacos" con evidencia real. |
| **CPIC** tablas de definición de alelos + guías | Alelos estrella CYP2C19/2D6/2C9/SLCO1B1/DPYD/TPMT/VKORC1… desde tag-SNPs del array | < 10 MB | Sólo alelos detectables en array (se marca lo no detectable). |
| **DGIdb** interactions.tsv + categorías | gen → fármacos/compuestos (aprobados e **investigacionales**) | ~30 MB | Semilla para "amplificar/modular" un gen. |
| **Open Targets Platform** (`molecule`, `knownDrugsAggregated`, `targets` — parquet) | Moléculas en fase 0-4 por diana, mecanismo de acción, descripción funcional del gen | ~300-500 MB | Mejor fuente abierta de compuestos experimentales por diana. |
| **dbSNP / gnomAD** | NO se descargan (decenas de GB). Frecuencias alélicas de los ~100 rsIDs priorizados se piden online a Ensembl REST (batch). | 0 | |

Total estimado: **~1.3-1.7 GB descargado**, DuckDB compactado ~0.7-1 GB. `fetch-data` mide tamaños reales, es idempotente/resumible, y aborta si excede 2 GB (indica qué omitir: p.ej. `targets` de Open Targets).

---

## 3. Análisis genético (fase 2: sólo datos locales)

### 3.1 Parseo
- Detección automática de formato. Normaliza a `(rsid, chr, pos, genotipo)`, build GRCh37. Descarta no-calls (`--`), maneja indels (`DI`, `II`, `DD`) y haploides (X en varones, Y, MT).
- **Hebra/orientación**: 23andMe reporta en hebra + de GRCh37. SNPedia indica orientación por SNP (plus/minus) → complementar genotipo cuando toque. GWAS Catalog da alelo de riesgo en hebra + (ambigüedades A/T, C/G se marcan "orientación incierta"). Punto crítico de corrección: test dedicado.

### 3.2 Cruce sistemático (exhaustivo)
Join de TODOS los SNPs del archivo contra SNPedia, ClinVar, GWAS, PharmGKB → tabla de hallazgos crudos (fuente, gen, efecto, medida de efecto, nivel de evidencia, genotipo del usuario, orientación resuelta).

### 3.3 Panel curado (profundidad donde el cruce automático es flojo)
YAML en el repo con ~250-400 variantes/haplotipos accionables con **tamaño de efecto y mecanismo** explícitos, por módulo funcional:

- **Metilación / un-carbono**: MTHFR C677T + A1298C (compuesto), MTRR A66G, MTR A2756G, MTHFD1, CBS, BHMT, SHMT1, PEMT (colina), FOLH1, TCN2, FUT2 (B12 sérica), SLC19A1.
- **Neurotransmisión / cognición / productividad**: COMT V158M (+ rs4633/rs4818 haplotipo), MAOA/MAOB, DRD2 Taq1A, DRD4, DAT1/SLC6A3, BDNF V66M, KIBRA, CHRNA4, SNAP25, TPH2, GAD1, GABRA2, ADRA2B, NRG1, CACNA1C, OXTR (5-HTTLPR: no en array, se indica).
- **Sueño / cafeína / cronotipo**: CYP1A2 (-163C>A), ADORA2A, PER2/PER3, CLOCK, RGS16, AANAT, MTNR1B.
- **Energía / mitocondria / estrés oxidativo**: PPARGC1A, NFE2L2 (NRF2), SOD2 A16V, GPX1, CAT, NQO1, GSTP1 (GSTM1/GSTT1 null: no detectable en array), UCP1/2/3, PPARA, PPARD, AMPD1, ACTN3, ACE I/D (tag).
- **Vitaminas / minerales**: VDR (Taq/Bsm/Fok), GC, CYP2R1, CYP24A1 (vit D); BCMO1 (β-caroteno→retinol); SLC23A1/2 (vit C); SLC30A8, SLC39A (zinc); HFE C282Y/H63D, TMPRSS6, TF (hierro); NBPF3/ALPL (B6); APOA5/CETP/LPL (lípidos); FADS1/FADS2 (conversión omega-3/6); TTPA (vit E); GGCX/VKORC1 (vit K); TRPM6/TRPM7/SHROOM3 (magnesio); SELENOP/GPX (selenio); PEMT/CHDH (colina); NAT2 (acetilación).
- **Inflamación / inmunidad / intestino**: IL6, IL1B, IL10, TNF, CRP, FUT2 (secretor, microbioma), HLA-DQ2/8 tag-SNPs (celiaquía), NOD2, ATG16L1, MCM6 (lactasa), AOC1/DAO y HNMT (histamina), ALDH2/ADH1B (alcohol).
- **Cardiometabólico / longevidad**: APOE (rs429358+rs7412 → ε2/ε3/ε4), LPA, PCSK9, LDLR, APOB, 9p21, TCF7L2, FTO, MC4R, PPARG, KLOTHO KL-VS, FOXO3, SIRT1/6, TERT, CDKN2B-AS, IGF1R, GHR, SH2B3, ADIPOQ, TP53 P72R, MTHFR↔homocisteína.
- **Farmacogenómica (CPIC/PharmGKB)**: CYP2C19, CYP2D6 (tags), CYP2C9, CYP3A4/5, CYP1A2, CYP2B6, SLCO1B1, VKORC1, DPYD, TPMT/NUDT15, UGT1A1, OPRM1, COMT (analgesia), G6PD, F5 Leiden, F2, SERPINA1, BCHE, RYR1/CACNA1S, ABCB1, ADRB1/2 (HLA-B: no detectable → se indica).
- **Rendimiento / composición corporal**: ACTN3 R577X, ACE, PPARGC1A, AMPD1, COL1A1/COL5A1 (tendón), IL6, MSTN, VDR.

Cada entrada: rsid(s), alelos de efecto con orientación, genotipos → fenotipo, magnitud de efecto (cuantitativa cuando exista: % actividad enzimática, OR, β), evidencia (meta-análisis / RCT / GWAS replicado / mecanístico), y **"palancas semilla"** (p.ej. MTHFR TT → L-metilfolato, riboflavina baja homocisteína ~40 % en TT; COMT Met/Met → cuidado con donantes de metilo/estimulantes; CYP1A2 lento → cafeína ↓; FADS1 minor → EPA/DHA preformado, no ALA; BCMO1 → retinol preformado; VDR/GC → D3+K2 dosis mayor y medir 25-OH-D; HFE → evitar hierro y vit C en exceso; APOE4 → DHA en fosfolípidos, vigilar Lp(a)…). Las semillas son punto de partida: el agente debe confirmarlas, matizarlas y **ampliarlas con lo novedoso**.

### 3.4 Haplotipos / combinaciones
Reglas YAML simples (no motor genérico) para APOE, CYP2C19 *2/*3/*17, CYP2D6 (tags detectables, con incertidumbre explícita), MTHFR compuesto heterocigoto, HFE compuesto, ACTN3×ACE, COMT×MAOA×DRD2 ("perfil dopaminérgico"), FADS1×FADS2, VDR×GC×CYP2R1 ("perfil vit D").

### 3.5 Limitaciones (se imprimen en el informe)
Arrays tienen tasa alta de falsos positivos en variantes raras (ClinVar patogénicas → "confirmar con secuenciación clínica"); sin imputación; sin CNV (GSTM1/GSTT1 null); sin HLA fiable; repeat-polymorphisms no detectables; orientación ambigua en SNPs A/T-C/G.

---

## 4. Significancia: umbrales explícitos (score.py)

**Hallazgo genético se reporta si cumple alguno:**
- SNPedia magnitude ≥ 2.5 para el genotipo del usuario (≥ 3 = destacado).
- ClinVar Pathogenic/Likely pathogenic con revisión ≥ 2 estrellas.
- GWAS: p ≤ 5e-8 **y** OR ≥ 1.3 (o ≤ 0.77) o |β| ≥ 0.1 SD, en rasgo de categoría objetivo; por rasgo se toma el estudio de mayor N.
- PharmGKB nivel 1A-2B / CPIC A-B.
- Panel curado: efecto funcional cuantificado (≥ 30 % cambio de actividad enzimática, fenotipo metabolizador alterado, o haplotipo de riesgo definido).

Score compuesto = f(tamaño de efecto, nivel de evidencia, accionabilidad, relevancia a objetivos [salud/energía/cognición-productividad/longevidad]). Ordena y **corta** (top ~30-50 hallazgos en ~8-14 módulos funcionales).

**Intervención se reporta si:**
- Evidencia humana: meta-análisis/RCT con efecto clínicamente relevante (cambio de biomarcador u outcome relevante, SMD ≥ ~0.3, o guía CPIC), **o**
- Compuesto novedoso/experimental: mecanismo directo sobre el gen/ruta afectada + al menos datos humanos preliminares o animales robustos, **etiquetado Experimental**, con estatus legal, riesgos y dosis de la literatura.
- Y **vinculada explícitamente** a genotipos del usuario (no "es bueno en general").

Grado de evidencia A-E impreso en cada intervención. Al agente se le exige excluir lo de efecto marginal ("estadísticamente significativo pero irrelevante").

---

## 5. Investigación con agente ACP (fase 3: online)

### 5.1 Mecánica ACP
- `acp_client.py`: `spawn_agent_process(client, *cmd)` → `initialize` → `new_session(cwd=tmp, mcp_servers=[])` → `prompt(...)`. `session_update` acumula texto y muestra progreso (spinner con último tool-call). `request_permission`: **permite** lecturas y herramientas web (WebSearch/WebFetch), **deniega** escritura/bash/edición. Timeout por tarea (15 min) + 1 reintento. Varias sesiones concurrentes sobre la misma conexión (3 por defecto; degrada a secuencial si el agente no lo soporta).
- El agente devuelve JSON en bloque ```json con esquema pydantic (`interventions[]`: name, class [suplemento/péptido/fármaco Rx/experimental/estilo de vida], target_genes[], rsids[], mechanism, expected_effect, evidence_grade, key_studies[] (PMID/DOI/NCT), dose, risks_interactions, novelty_flag, legal_status; `labs_to_order[]`). Parse → reintento con corrección si falla.
- Sólo viajan al agente: rsids priorizados + genotipos + extractos de evidencia local. Nunca el archivo.

### 5.2 Tareas
1. **Por módulo funcional** (≈ 8-14 prompts en paralelo): contexto = genotipos del módulo, fenotipo inferido ("actividad MTHFR ~30 %", "metabolizador lento CYP1A2"), extractos SNPedia/PharmGKB/DGIdb/Open Targets (fármacos conocidos e investigacionales para las dianas del módulo), palancas semilla. Instrucciones: buscar exhaustivamente intervenciones con efecto significativo para ESTE perfil; cubrir obligatoriamente suplementos, péptidos, fármacos (Rx incluidos, etiquetados), compuestos experimentales/novedosos (literatura 2023-2026, preprints bioRxiv/medRxiv, ClinicalTrials.gov); mecanismo → gen, dosis, evidencia, riesgos; descartar lo marginal; "qué medir en analítica" para el módulo (homocisteína, 25-OH-D, B12/MMA, ferritina, ApoB, Lp(a), hs-CRP, índice omega-3…).
2. **Auditoría de evidencia (determinista, sin LLM)**: para cada compuesto propuesto, `evidence.py` consulta Europe PMC (nº de RCT/meta-análisis humanos, últimos 5 años, títulos top) y ClinicalTrials.gov (ensayos activos). Se adjunta y se rebaja/etiqueta "experimental" si no hay datos humanos. Ancla la salida del agente y reduce alucinación.
3. **Síntesis final** (1 prompt con todo): stack integrado priorizado, **conflictos e interacciones entre módulos** (COMT lento × donantes de metilo; HFE × vit C alta; APOE4 × ciertas grasas; CYP2C19 PM × PPIs/clopidogrel; G6PD × vit C IV/azul de metileno; SLCO1B1 × estatinas; CYP1A2 lento × cafeína+estimulantes; ADRA2B/COMT × nootrópicos dopaminérgicos), orden de introducción, qué confirmar con analítica/secuenciación, y hallazgos "inesperados" del cruce automático fuera del panel.

### 5.3 Sin agente
Si el agente no arranca (o `--offline`), se imprime el análisis local completo (hallazgos, módulos, palancas semilla) con aviso. Misma salida sin la sección de investigación; no es una feature extra.

---

## 6. Informe (consola + `.md`)

1. Cabecera: formato detectado, nº SNPs, nº cruzados, build, limitaciones.
2. **Resumen ejecutivo** (10-15 líneas: lo que más importa y el stack central).
3. Por módulo funcional: tabla de genotipos clave (rsid/gen/genotipo/fenotipo/evidencia), interpretación, **intervenciones** ordenadas por impacto × evidencia con grado A-E, dosis, mecanismo, riesgos, flag Experimental/Rx, refs (PMID/DOI/NCT).
4. Hallazgos significativos fuera del panel (SNPedia mag ≥ 3, ClinVar, GWAS fuertes) con nota de confirmación.
5. Farmacogenómica (fenotipos metabolizadores y fármacos afectados).
6. Interacciones/conflictos del stack + orden de introducción + analíticas a pedir.
7. Fuentes y versiones de bases usadas; disclaimer (no es consejo médico; compuestos experimentales: estatus legal y seguridad).

---

## 7. Hitos

| Hito | Contenido | Criterio de hecho |
|---|---|---|
| M0 | Scaffold, `fetch-data` con presupuesto 2 GB, DuckDB, descarga resumible | `genostack fetch-data` termina < 2 GB y reporta tamaños |
| M1 | Parsers + normalización de hebra + cruce contra todas las bases | Fixture sintético 23andMe → tabla de hallazgos crudos correcta (tests de orientación) |
| M2 | Panel YAML + haplotipos + scoring/filtro + módulos | Informe offline completo y legible |
| M3 | Cliente ACP (permisos, streaming, concurrencia, JSON) + prompts por módulo + auditoría Europe PMC/CT.gov + síntesis | `analyze` end-to-end con `claude-agent-acp`; test con agente mock |
| M4 | Pulido del informe, disclaimer, `.md`, README | Revisión final con un genoma real |

Tamaño estimado: ~2.500-3.500 líneas Python + ~1.500 líneas YAML de panel.

---

## 8. Supuestos (corrígelos si no encajan)

- "sub por ACP" = subagente/agente externo vía Agent Client Protocol; cliente = nuestro programa, agente por defecto = Claude Code (`claude-agent-acp`), intercambiable con `--agent`.
- Idioma del informe: español.
- Objetivos ponderados por defecto: salud general, energía, cognición/productividad, longevidad.
- Se incluyen fármacos con receta y compuestos experimentales/péptidos, siempre etiquetados con estatus y riesgos.
- Salida por consola; además se guarda `.md` al lado del input porque el texto es largo.

---

## 9. Estado de implementación (2026-08-23)

Hecho y probado:
- `genostack fetch-data`: descarga resumible (ClinVar 442 MB, GWAS Catalog 73 MB, DGIdb 12 MB, Open Targets ~110 MB, CPIC API, SNPedia 216k páginas vía API en 2 fases paralelas) → `data/genostack.duckdb` (186 MB; 0.82 GB descargados, presupuesto 2 GB respetado). PharmGKB directo ya no existe (`api.pharmgkb.org` desapareció): se usa el dataset `pharmacogenomics` de Open Targets (derivado de PharmGKB/ClinPGx) + CPIC.
- Parsers (23andMe/Ancestry/MyHeritage/FTDNA/VCF), alias de ids internos 23andMe (`iNNNNNN` → rsid vía SNPedia), cruce SNPedia (con orientación de hebra) / ClinVar / GWAS (filtros de unidades, ruido pQTL/NMR, alelo mayor) / PGx (genotipos de referencia excluidos), panel curado de **207 variantes + 13 reglas de haplotipo** en 9 módulos (alelos verificados contra Ensembl/VEP/SNPedia por los subagentes), scoring y agrupación por módulos, informe Markdown en consola + `.md`.
- Cliente ACP (`acp_client.py`): spawn del agente, permisos (sólo lectura/búsqueda/web), sesiones concurrentes, timeouts, extracción de JSON, logs de prompts/respuestas. Probado end-to-end con un agente mock (`tests/mock_agent.py`) incluyendo auditoría de evidencia real en Europe PMC / ClinicalTrials.gov y síntesis.
- Tests: `python -m pytest -q` (10 tests).

Pendiente del lado del usuario:
- El adaptador `claude-agent-acp` no pudo autenticarse en esta máquina ("OAuth session expired"): hay que hacer login en Claude Code desde una terminal (`claude` → `/login`) o exportar `ANTHROPIC_API_KEY`; después `genostack analyze <raw.txt>` ejecuta la investigación real (≈9 módulos en paralelo + síntesis; 10-30 min según el agente).
