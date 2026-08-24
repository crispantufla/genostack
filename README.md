# genostack

Analizador local de raw reports de ADN (23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA, VCF simple) que:

1. cruza **todos** los SNPs del archivo contra bases locales descargadas de antemano (≤ 2 GB): SNPedia, ClinVar, GWAS Catalog, Open Targets (incl. farmacogenómica PharmGKB/ClinPGx), DGIdb y CPIC, más un **panel curado** de variantes/haplotipos accionables (metilación, neurotransmisión, sueño/cafeína, mitocondria, vitaminas/minerales, inflamación/intestino, cardiometabólico/longevidad, rendimiento, farmacogenómica);
2. filtra por **significancia** (umbrales explícitos: SNPedia magnitud ≥ 2.5, ClinVar ≥ 2★, GWAS p ≤ 5e-8 con OR ≥ 1.3 o β relevante, PharmGKB 1A-2B, panel con efecto funcional cuantificado);
3. delega a un **agente por ACP** (Agent Client Protocol; por defecto Claude Code vía `claude-agent-acp`) la investigación exhaustiva de intervenciones —suplementos, péptidos, fármacos (Rx), compuestos experimentales/novedosos 2023-2026— vinculadas a los genotipos concretos, con una **auditoría determinista** de evidencia en Europe PMC / ClinicalTrials.gov y una síntesis final (stack, conflictos, orden, analíticas);
4. imprime el informe en consola (y guarda una copia `.md` junto al input).

## Uso

```bash
uv sync                                     # crea el entorno e instala (https://docs.astral.sh/uv/)

uv run genostack fetch-data                 # una vez; reanudable; ~0.9 GB descargados
uv run genostack analyze genome_23andme.txt # informe completo (requiere agente ACP autenticado)
uv run genostack analyze genome_23andme.txt --offline   # sólo análisis local (sin agente ni red)
```

Opciones de `analyze`: `--agent "<cmd>"` (otro agente ACP: `gemini --experimental-acp`, `codex-acp`…), `--concurrency N`, `--goals health,energy,cognition,longevity,performance,pharmaco`, `--lang es|en`, `--no-save`.

### Agente ACP

Por defecto se lanza `npx -y @agentclientprotocol/claude-agent-acp` (necesita Node ≥ 18 y una sesión de Claude Code iniciada: ejecuta `claude` en una terminal y haz `/login`, o exporta `ANTHROPIC_API_KEY`). El programa actúa como *cliente* ACP: concede al agente sólo herramientas de lectura/búsqueda/web y deniega escritura/ejecución. Al agente sólo viajan los rsIDs ya priorizados con su genotipo y extractos de evidencia local; el archivo crudo nunca sale de la máquina. Los prompts y respuestas se guardan en un directorio temporal (`genostack-agent-*`) indicado al final del informe.

## Datos locales

`fetch-data` descarga a `./data/raw` (o `$GENOSTACK_DATA`) y construye `data/genostack.duckdb`:

| Fuente | Tablas | Notas |
|---|---|---|
| SNPedia (API MediaWiki, CC-BY-NC-SA) | `snpedia_snp`, `snpedia_geno` | magnitud/repute por genotipo, orientación de hebra |
| ClinVar `variant_summary` | `clinvar` | P/LP, factor de riesgo, respuesta a fármaco, protectora; estrellas de revisión |
| GWAS Catalog (asociaciones completas) | `gwas` | p ≤ 5e-8, alelo de riesgo, OR/β, N |
| Open Targets Platform | `ot_molecule`, `ot_moa`, `ot_clinical_target`, `ot_target`, `ot_pgx` | fármacos por diana (fase 0-4), mecanismo, función génica, PGx (PharmGKB) |
| DGIdb | `dgidb` | interacciones gen-fármaco (aprobados e investigacionales) |
| CPIC API | `cpic_pair`, `cpic_allele` | pares gen-fármaco nivel A/B, alelos estrella |

Presupuesto máximo 2 GB (se comprueba). El panel curado vive en `genostack/panel/*.yaml` (esquema en `panel/SCHEMA.md`; alelos siempre en hebra + GRCh37).

## Limitaciones

Arrays de consumo: sin imputación, alta tasa de falsos positivos en variantes raras (confirmar ClinVar patogénicas con secuenciación clínica), sin CNV (GSTM1/GSTT1), sin tipado HLA real, SNPs A/T-C/G con orientación incierta marcados. No es consejo médico.

## Tests

```bash
uv run pytest -q
```
