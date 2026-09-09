"""Predicción del efecto molecular de una variante con AlphaGenome (DeepMind).

Qué añade y qué NO añade
------------------------
El resto del pipeline trabaja con evidencia **empírica**: asociaciones GWAS medidas en poblaciones,
anotaciones clínicas de ClinVar, guías CPIC, ensayos. AlphaGenome no observa nada: **predice** desde
la secuencia qué le hace una variante a la maquinaria de la célula (expresión génica, splicing,
accesibilidad de cromatina) en cientos de tipos celulares.

Por eso aquí sólo se usa para una cosa concreta: **poner mecanismo donde no lo hay**. Un hallazgo GWAS
intergénico dice "esto se asocia a tal rasgo" sin decir sobre qué gen actúa ni cómo; AlphaGenome puede
sugerir "reduce la expresión de tal gen en hígado". Para las variantes del panel curado (MTHFR, APOE,
PNPLA3...) el mecanismo lleva descrito décadas y esto no aporta.

Reglas que se respetan en todo el programa
------------------------------------------
1. Sus resultados se etiquetan SIEMPRE como predicción computacional, nunca como observación.
2. **No pueden subir el grado de evidencia** de ninguna intervención ni el score de ningún hallazgo.
3. DeepMind declara explícitamente que AlphaGenome «no ha sido validado ni está aprobado para ningún
   uso clínico», y que sirve como parte de una cadena de evidencia, nunca como evidencia suficiente.

Requisitos
----------
* `pip install alphagenome` (o `uv sync --extra alphagenome`).
* Clave de API gratuita para uso no comercial: https://deepmind.google.com/science/alphagenome
  Se lee de la variable de entorno `ALPHAGENOME_API_KEY`.
* Las coordenadas del modelo son GRCh38; los raw reports son GRCh37, así que cada rsID se resuelve
  antes contra Ensembl (que devuelve GRCh38 y el alelo de referencia).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import httpx

ENSEMBL = "https://rest.ensembl.org/variation/homo_sapiens"
INTERVAL_WIDTH = 131_072          # ventana centrada en la variante (SEQUENCE_LENGTH_128KB)
MAX_SCORERS = 20                  # tope del servicio por petición

# Consecuencias donde el mecanismo YA se conoce: no gastamos predicciones ahí.
_WELL_UNDERSTOOD = {
    "missense_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "inframe_deletion", "inframe_insertion", "synonymous_variant",
}


@dataclass
class Prediction:
    rsid: str
    chrom: str
    position: int
    ref: str
    alt: str
    consequence: str = ""
    effects: list[dict] = field(default_factory=list)   # efectos más marcados, ya ordenados
    error: str | None = None

    def summary(self, n: int = 3) -> str:
        if self.error:
            return f"sin predicción ({self.error})"
        if not self.effects:
            return "sin efecto molecular destacable"
        out = []
        for e in self.effects[:n]:
            gene = e.get("gene") or "?"
            what = e.get("output") or "señal"
            tissue = f" en {e['tissue']}" if e.get("tissue") else ""
            direction = "↑" if (e.get("score") or 0) > 0 else "↓"
            out.append(f"{direction} {what} de {gene}{tissue} (q={e.get('quantile', 0):.2f})")
        return "; ".join(out)


# ----------------------------------------------------------------------------- coordenadas
def resolve_grch38(rsids: list[str], timeout: float = 90) -> dict[str, dict]:
    """rsID -> {chrom, pos, ref, alts, consequence} en GRCh38, vía Ensembl REST (lotes de 200)."""
    rsids = [r for r in dict.fromkeys(rsids) if r.startswith("rs")]
    out: dict[str, dict] = {}
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    with httpx.Client(timeout=timeout, headers=headers) as client:
        for i in range(0, len(rsids), 200):
            chunk = rsids[i:i + 200]
            try:
                r = client.post(ENSEMBL, json={"ids": chunk})
                if r.status_code != 200:
                    continue
                for rsid, v in r.json().items():
                    maps = v.get("mappings") or []
                    if not maps:
                        continue
                    m = maps[0]
                    alleles = (m.get("allele_string") or "").split("/")
                    if len(alleles) < 2 or any(len(a) != 1 for a in alleles):
                        continue          # indels y multialélicos complejos: fuera
                    out[rsid] = {
                        "chrom": str(m.get("seq_region_name")),
                        "pos": int(m.get("start")),
                        "ref": alleles[0].upper(),
                        "alts": [a.upper() for a in alleles[1:]],
                        "consequence": v.get("most_severe_consequence") or "",
                    }
            except (httpx.HTTPError, ValueError, KeyError):
                continue
    return out


def pick_alt(user_genotype: str, ref: str, alts: list[str]) -> str | None:
    """El alelo del usuario que NO es el de referencia; None si es homocigoto de referencia."""
    carried = {c for c in user_genotype.upper() if c in "ACGT"} - {ref}
    for a in alts:
        if a in carried:
            return a
    return None


def worth_predicting(consequence: str) -> bool:
    """Sólo tiene sentido predecir donde el mecanismo no es evidente por la consecuencia."""
    return consequence not in _WELL_UNDERSTOOD


# ----------------------------------------------------------------------------- predicción
def predict(variants: list[tuple[str, str]], api_key: str | None = None,
            progress=None, limit: int = 60) -> dict[str, Prediction]:
    """Puntúa (rsid, genotipo_del_usuario) y devuelve los efectos moleculares más marcados.

    `limit` acota cuántas variantes se envían: cada una es una petición al servicio y sólo interesan
    las priorizadas, no las 600.000 del array.
    """
    key = api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        raise RuntimeError(
            "Falta la clave de AlphaGenome. Consíguela gratis para uso no comercial en "
            "https://deepmind.google.com/science/alphagenome y expórtala como ALPHAGENOME_API_KEY.")
    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client, variant_scorers
    except ImportError as e:
        raise RuntimeError("Falta el paquete `alphagenome` (uv sync --extra alphagenome)") from e

    coords = resolve_grch38([r for r, _ in variants])
    out: dict[str, Prediction] = {}
    todo: list = []
    for rsid, gt in variants:
        c = coords.get(rsid)
        if not c:
            out[rsid] = Prediction(rsid, "", 0, "", "", error="sin coordenadas en Ensembl")
            continue
        alt = pick_alt(gt, c["ref"], c["alts"])
        if alt is None:
            continue                            # homocigoto de referencia: nada que predecir
        p = Prediction(rsid, c["chrom"], c["pos"], c["ref"], alt, c["consequence"])
        if not worth_predicting(c["consequence"]):
            p.error = f"mecanismo ya conocido ({c['consequence']})"
            out[rsid] = p
            continue
        out[rsid] = p
        todo.append(p)

    todo = todo[:limit]
    if not todo:
        return out

    client = dna_client.create(key)
    scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())[:MAX_SCORERS]
    for n, p in enumerate(todo, 1):
        try:
            variant = genome.Variant(chromosome=f"chr{p.chrom}", position=p.position,
                                     reference_bases=p.ref, alternate_bases=p.alt, name=p.rsid)
            half = INTERVAL_WIDTH // 2
            interval = genome.Interval(chromosome=f"chr{p.chrom}",
                                       start=max(0, p.position - half), end=p.position + half)
            scores = client.score_variant(interval=interval, variant=variant, variant_scorers=scorers)
            p.effects = _top_effects(variant_scorers.tidy_scores(scores))
        except Exception as e:                                        # noqa: BLE001
            p.error = f"{type(e).__name__}: {e!s:.120}"
        if progress:
            progress(n, len(todo), p)
    return out


def _top_effects(df, n: int = 6) -> list[dict]:
    """Se queda con los efectos más extremos del DataFrame que devuelve `tidy_scores`.

    Se leen las columnas de forma defensiva: el esquema exacto depende de la versión del SDK, y lo
    único garantizado es que hay una fila por score con alguna medida de magnitud.
    """
    if df is None or len(df) == 0:
        return []
    cols = {c.lower(): c for c in df.columns}
    score_col = next((cols[c] for c in ("quantile_score", "raw_score", "score") if c in cols), None)
    if score_col is None:
        return []
    gene_col = next((cols[c] for c in ("gene_name", "gene_symbol", "gene_id", "gene") if c in cols), None)
    out_col = next((cols[c] for c in ("output_type", "output", "assay") if c in cols), None)
    tis_col = next((cols[c] for c in ("biosample_name", "ontology_curie", "biosample", "track_name") if c in cols), None)
    d = df.assign(_abs=df[score_col].abs()).sort_values("_abs", ascending=False).head(n)
    effects = []
    for _, row in d.iterrows():
        effects.append({
            "gene": str(row[gene_col]) if gene_col else None,
            "output": str(row[out_col]) if out_col else None,
            "tissue": str(row[tis_col]) if tis_col else None,
            "score": float(row[score_col]),
            "quantile": float(row[cols["quantile_score"]]) if "quantile_score" in cols else float(row[score_col]),
        })
    return effects
