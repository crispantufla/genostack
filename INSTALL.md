# genostack — Instalación y uso

Analiza tu raw report de ADN (23andMe, AncestryDNA, MyHeritage, FTDNA) en tu propio PC y genera un informe profundo: hallazgos genéticos significativos + intervenciones investigadas por un agente de IA (suplementos, péptidos, fármacos, compuestos experimentales) vinculadas a tu genotipo. **Tu genoma nunca sale de tu ordenador**; al agente solo se le envían los ~30-80 marcadores ya priorizados.

> ⚠ Esto no es consejo médico. Los arrays de consumo tienen falsos positivos en variantes raras: cualquier hallazgo patogénico debe confirmarse con secuenciación clínica, y los fármacos/compuestos experimentales deben discutirse con un profesional.

## 1. Requisitos

- **Windows, macOS o Linux**
- **uv** — https://docs.astral.sh/uv/ (Windows: `winget install astral-sh.uv` · macOS/Linux: `curl -LsSf https://astral.sh/uv/install.sh | sh`). uv instala solo el Python que haga falta.
- **Node.js 18 o superior** — https://nodejs.org (solo para el agente de IA)
- ~1 GB de disco para las bases de datos locales
- Una cuenta de Claude (Pro/Max o clave de API de Anthropic) para la fase de investigación

## 2. Instalación

Descomprime el zip donde quieras (p.ej. `C:\genostack`) y en una terminal (PowerShell en Windows):

```powershell
cd C:\genostack
uv sync
```

Eso crea el entorno e instala todo. A partir de ahí, todos los comandos van con `uv run` (sin activar nada).

<details><summary>Sin uv (pip clásico)</summary>

```powershell
python -m venv .venv
.\.venv\Scripts\activate          # Linux/macOS: source .venv/bin/activate
pip install -e .
```
y usa `genostack ...` en lugar de `uv run genostack ...`.
</details>

## 3. Descargar las bases de datos locales (una sola vez)

```powershell
uv run genostack fetch-data
```

Descarga ~0.8 GB (SNPedia, ClinVar, GWAS Catalog, Open Targets, DGIdb, CPIC) y construye `data\genostack.duckdb`. Tarda 30-60 min (la mayor parte es SNPedia). **Es reanudable**: si se corta, vuelve a ejecutarlo y continúa donde iba.

Atajo: si alguien te pasa su carpeta `data\` ya construida, cópiala junto al programa y sáltate este paso.

## 4. Autenticar el agente de IA (una sola vez)

El agente (Claude Code) se ejecuta localmente y usa TU cuenta. Opción A (recomendada, usa tu suscripción):

```powershell
npm install -g @anthropic-ai/claude-code
claude auth login
```

Opción B (clave de API, se cobra por uso): antes de cada análisis ejecuta

```powershell
$env:ANTHROPIC_API_KEY = "sk-ant-..."     # Linux/macOS: export ANTHROPIC_API_KEY=sk-ant-...
```

## 5. Analizar tu genoma

Descarga tu "raw data" desde tu proveedor (23andMe: Cuenta → Browse Raw Data → Download; Ancestry: Configuración → Descargar datos de ADN). Sirve el `.txt`, `.csv`, `.zip` o `.gz` tal cual. Después:

```powershell
uv run genostack analyze "C:\Users\tu\Downloads\genome_XXXX.zip"
```

- Fase local (~20 s): cruce contra las bases y el panel curado.
- Fase de investigación (10-30 min): ~9 módulos en paralelo + auditoría de evidencia (Europe PMC / ClinicalTrials.gov) + síntesis final. Verás el progreso en pantalla.
- El informe se imprime en consola y se guarda como `<tu_archivo>.genostack.md`.

### Opciones útiles

| Opción | Qué hace |
|---|---|
| `--offline` | Solo análisis local, sin agente ni red (informe sin la sección de intervenciones) |
| `--goals health,energy,cognition,longevity,performance,pharmaco` | Pondera los objetivos (por defecto: health,energy,cognition,longevity) |
| `--concurrency 3` | Sesiones de investigación simultáneas (súbelo a 4-5 si tu plan lo aguanta) |
| `--lang en` | Informe en inglés |
| `--agent "gemini --experimental-acp"` | Usa otro agente ACP (Gemini CLI, codex-acp…) en lugar de Claude |
| `--no-save` | No guardar el `.md`, solo consola |

## 6. Problemas frecuentes

- **"El agente ACP no pudo autenticarse / OAuth session expired"** → ejecuta `claude auth login` (o define `ANTHROPIC_API_KEY`) y reintenta.
- **"No existe la base local"** → ejecuta `genostack fetch-data`.
- **`genostack` no se reconoce** → usa `uv run genostack ...` desde la carpeta del proyecto (o activa el entorno `.venv`).
- **npx/node no encontrado** → instala Node.js 18+ y reabre la terminal.
- **SNPedia va lento o da 502** → es normal; el proceso reintenta y es reanudable.
- Los prompts y respuestas del agente quedan en un directorio `genostack-agent-*` (la ruta sale al final del informe) por si quieres auditar qué se envió: solo rsIDs priorizados, nunca el archivo completo.

## 7. Verificar la instalación (opcional)

```powershell
uv run pytest -q             # 10 tests (uno se salta si aún no ejecutaste fetch-data)
uv run genostack analyze tests\fixture_23andme.txt --offline   # informe de prueba con un genoma sintético (requiere fetch-data)
```
