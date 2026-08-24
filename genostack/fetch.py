"""`genostack fetch-data`: download the offline knowledge bases (<= 2 GB) and build one DuckDB file.

Sources
  - SNPedia (MediaWiki API; SNP pages + genotype pages)      -> snpedia_snp, snpedia_geno
  - ClinVar variant_summary.txt.gz                            -> clinvar
  - GWAS Catalog associations (full, ontology-annotated)      -> gwas
  - Open Targets: drug_molecule, drug_mechanism_of_action,
                  clinical_target, pharmacogenomics, target   -> ot_molecule, ot_moa, ot_clinical_target, ot_pgx, ot_target
  - DGIdb interactions.tsv                                    -> dgidb
  - CPIC API (gene-drug pairs, alleles, drugs)                -> cpic_pair, cpic_allele
Everything is resumable: re-running only fetches what is missing.
"""
from __future__ import annotations

import json
import re
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import httpx
from rich.console import Console
from rich.progress import (BarColumn, DownloadColumn, Progress, TextColumn, TimeRemainingColumn,
                           TransferSpeedColumn)

from .config import DATA_DIR, DB_PATH, DOWNLOAD_BUDGET_BYTES, MANIFEST_PATH, RAW_DIR

console = Console()

OT_BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output"
OT_DATASETS = ["drug_molecule", "drug_mechanism_of_action", "clinical_target", "pharmacogenomics", "target"]

FILES = {
    "clinvar": ("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz", "clinvar_variant_summary.txt.gz"),
    "gwas": ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip", "gwas_associations_full.zip"),
    "dgidb": ("https://dgidb.org/data/latest/interactions.tsv", "dgidb_interactions.tsv"),
}
CPIC_API = "https://api.cpicpgx.org/v1"
SNPEDIA_API = "https://bots.snpedia.com/api.php"
UA = "genostack/0.1 (local personal genome analysis; contact: none)"


# ----------------------------------------------------------------------------- helpers
def _manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    return {"sources": {}}


def _save_manifest(m: dict) -> None:
    MANIFEST_PATH.write_text(json.dumps(m, indent=2, ensure_ascii=False), encoding="utf-8")


def _raw_bytes() -> int:
    return sum(p.stat().st_size for p in RAW_DIR.rglob("*") if p.is_file())


def _budget_check() -> None:
    used = _raw_bytes()
    if used > DOWNLOAD_BUDGET_BYTES:
        raise SystemExit(f"Presupuesto de descarga superado: {used/1e9:.2f} GB > 2 GB. Borra data/raw parcial o ajusta fuentes.")


def download(url: str, dest: Path, client: httpx.Client, label: str | None = None) -> Path:
    """Resumable download with progress bar. Skips if a complete file exists."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    headers = {}
    existing = part.stat().st_size if part.exists() else 0
    if existing:
        headers["Range"] = f"bytes={existing}-"
    with client.stream("GET", url, headers=headers, follow_redirects=True, timeout=120) as r:
        if r.status_code == 416:  # range not satisfiable -> already complete
            part.rename(dest)
            return dest
        r.raise_for_status()
        if r.status_code != 206:
            existing = 0
            mode = "wb"
        else:
            mode = "ab"
        total = int(r.headers.get("Content-Length", 0)) + existing or None
        with Progress(TextColumn("[bold blue]{task.fields[name]}"), BarColumn(), DownloadColumn(),
                      TransferSpeedColumn(), TimeRemainingColumn(), console=console, transient=True) as prog:
            task = prog.add_task("dl", total=total, completed=existing, name=label or dest.name)
            with open(part, mode) as f:
                for chunk in r.iter_bytes(1 << 20):
                    f.write(chunk)
                    prog.update(task, advance=len(chunk))
    part.rename(dest)
    return dest


def _http() -> httpx.Client:
    return httpx.Client(headers={"User-Agent": UA}, timeout=120)


# ----------------------------------------------------------------------------- SNPedia
def fetch_snpedia(client: httpx.Client, workers: int = 4) -> None:
    """Scrape SNP pages and genotype pages via the MediaWiki API into JSONL caches (resumable, parallel).

    Phase A lists all titles of the category (500 per request); phase B fetches page content in batches of 50
    titles with a small thread pool. Titles already present in the JSONL cache are skipped.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    out_dir = RAW_DIR / "snpedia"
    out_dir.mkdir(parents=True, exist_ok=True)
    for cat, fname in (("Category:Is_a_snp", "snps.jsonl"), ("Category:Is_a_genotype", "genotypes.jsonl")):
        jsonl = out_dir / fname
        state_f = out_dir / (fname + ".state")
        titles_f = out_dir / (fname + ".titles")
        if state_f.exists() and state_f.read_text().strip() == "DONE":
            continue
        # ---- phase A: titles
        if titles_f.exists() and titles_f.read_text(encoding="utf-8").endswith("\n#DONE\n"):
            titles = [t for t in titles_f.read_text(encoding="utf-8").splitlines() if t and not t.startswith("#")]
        else:
            titles = []
            cont: dict = {}
            total = _snpedia_cat_size(client, cat)
            console.print(f"[cyan]SNPedia[/] {cat}: listando {total:,} títulos…")
            with open(titles_f, "w", encoding="utf-8") as tf:
                while True:
                    d = _snpedia_get(client, {"action": "query", "list": "categorymembers", "cmtitle": cat, "cmlimit": "500",
                                              "cmprop": "title", "format": "json", "formatversion": "2", **cont})
                    batch = [m["title"] for m in d.get("query", {}).get("categorymembers", [])]
                    titles.extend(batch)
                    tf.write("\n".join(batch) + "\n")
                    if "continue" not in d:
                        tf.write("#DONE\n")
                        break
                    cont = dict(d["continue"])
        # ---- phase B: content
        have: set[str] = set()
        if jsonl.exists():
            with open(jsonl, encoding="utf-8") as fh:
                for line in fh:
                    try:
                        have.add(json.loads(line)["t"])
                    except Exception:  # noqa: BLE001
                        continue
        todo = [t for t in titles if t not in have]
        console.print(f"[cyan]SNPedia[/] {cat}: {len(have):,}/{len(titles):,} páginas en caché; faltan {len(todo):,}")
        batches = [todo[i:i + 50] for i in range(0, len(todo), 50)]

        def fetch_batch(batch: list[str]) -> list[dict]:
            d = _snpedia_get(client, {"action": "query", "titles": "|".join(batch), "prop": "revisions", "rvprop": "content",
                                      "rvslots": "main", "format": "json", "formatversion": "2"})
            out = []
            for pg in d.get("query", {}).get("pages", []):
                revs = pg.get("revisions")
                if not revs:
                    continue
                out.append({"t": pg["title"], "c": revs[0].get("slots", {}).get("main", {}).get("content", "")})
            return out

        with Progress(TextColumn("[bold blue]{task.fields[name]}"), BarColumn(), TextColumn("{task.completed}/{task.total}"),
                      TimeRemainingColumn(), console=console, transient=True) as prog:
            task = prog.add_task("snpedia", total=len(titles), completed=len(have), name=cat.split(":")[1])
            with open(jsonl, "a", encoding="utf-8") as out, ThreadPoolExecutor(max_workers=workers) as ex:
                futs = {ex.submit(fetch_batch, b): b for b in batches}
                for fut in as_completed(futs):
                    recs = fut.result()
                    for r in recs:
                        out.write(json.dumps(r, ensure_ascii=False) + "\n")
                    out.flush()
                    prog.update(task, advance=len(futs[fut]))
        state_f.write_text("DONE")


def _snpedia_cat_size(client: httpx.Client, cat: str) -> int:
    d = _snpedia_get(client, {"action": "query", "prop": "categoryinfo", "titles": cat, "format": "json", "formatversion": "2"})
    return d["query"]["pages"][0]["categoryinfo"]["pages"]


def _snpedia_get(client: httpx.Client, params: dict) -> dict:
    for attempt in range(8):
        try:
            r = client.get(SNPEDIA_API, params=params, timeout=90)
            r.raise_for_status()
            return r.json()
        except (httpx.HTTPError, ValueError) as e:
            wait = min(60, 2 ** attempt)
            console.print(f"[yellow]SNPedia: {e!s:.80} — reintento en {wait}s[/]")
            time.sleep(wait)
    raise SystemExit("SNPedia no responde; vuelve a ejecutar fetch-data más tarde (es reanudable).")


_TEMPLATE_RE = re.compile(r"\{\{\s*(Rsnum|Genotype|23andMe SNP)\s*(.*?)\}\}", re.S | re.I)


def parse_template(content: str) -> tuple[str, dict[str, str]]:
    m = _TEMPLATE_RE.search(content)
    if not m:
        return "", {}
    body = m.group(2)
    fields: dict[str, str] = {}
    for part in re.split(r"\n\s*\|", "\n" + body):
        part = part.strip().lstrip("|").strip()
        if "=" in part:
            k, v = part.split("=", 1)
            fields[k.strip().lower()] = v.strip()
    return m.group(1).lower(), fields


_WIKI_STRIP = [
    (re.compile(r"\{\{[^{}]*\}\}", re.S), " "),          # templates
    (re.compile(r"\[\[(?:[^|\]]*\|)?([^\]]*)\]\]"), r"\1"),  # links
    (re.compile(r"\[https?://\S+\s?([^\]]*)\]"), r"\1"),
    (re.compile(r"<ref[^>]*/>|<ref[^>]*>.*?</ref>", re.S), " "),
    (re.compile(r"<[^>]+>"), " "),
    (re.compile(r"'{2,}"), ""),
    (re.compile(r"={2,}[^=]+={2,}"), " "),
    (re.compile(r"[ \t]+"), " "),
    (re.compile(r"\n{2,}"), "\n"),
]


def strip_wiki(text: str, limit: int = 1500) -> str:
    for rx, rep in _WIKI_STRIP:
        text = rx.sub(rep, text)
    text = text.strip()
    return text[:limit]


def build_snpedia_tables(con: duckdb.DuckDBPyConnection) -> None:
    snp_rows, geno_rows = [], []
    snps_f = RAW_DIR / "snpedia" / "snps.jsonl"
    geno_f = RAW_DIR / "snpedia" / "genotypes.jsonl"
    alias_rows = []
    for line in open(snps_f, encoding="utf-8"):
        rec = json.loads(line)
        title = rec["t"]
        if title.lower().startswith("i") and title[1:].isdigit():
            # 23andMe internal ids (iNNNNNN) -> rsid, from the {{23andMe SNP|rsid=...}} template
            _, f = parse_template(rec["c"])
            rs = (f.get("rsid") or "").strip().lower()
            if rs.isdigit():
                rs = "rs" + rs
            if not rs.startswith("rs"):
                m = re.search(r"\[\[(rs\d+)\]\]", rec["c"], re.I)
                rs = m.group(1).lower() if m else ""
            if rs.startswith("rs"):
                alias_rows.append((title.lower(), rs))
            continue
        if not title.lower().startswith("rs"):
            continue
        kind, f = parse_template(rec["c"])
        orient = (f.get("stabilizedorientation") or f.get("orientation") or "").lower()
        text = strip_wiki(rec["c"])
        snp_rows.append((title.lower(), f.get("gene") or f.get("gene_s") or None, f.get("chromosome"),
                         _int(f.get("position")), orient or None, f.get("summary") or None,
                         _float(f.get("gmaf")), text or None))
    for line in open(geno_f, encoding="utf-8"):
        rec = json.loads(line)
        kind, f = parse_template(rec["c"])
        if kind != "genotype":
            continue
        rsid = (f.get("rsid") or "").strip().lower()
        if rsid.isdigit():
            rsid = "rs" + rsid
        m = re.match(r"(rs\d+)\(([^;]+);([^)]+)\)", rec["t"].lower().replace(" ", ""))
        if m and (not rsid.startswith("rs") or rsid != m.group(1)):
            rsid = m.group(1)  # the page title is authoritative
        if not rsid.startswith("rs"):
            continue
        a1 = (f.get("allele1") or (m.group(2).upper() if m else "")).upper()
        a2 = (f.get("allele2") or (m.group(3).upper() if m else "")).upper()
        geno_rows.append((rsid, a1, a2, _float(f.get("magnitude")), (f.get("repute") or "").strip().lower() or None,
                          f.get("summary") or None))
    import pyarrow as pa

    def _arrow(rows, cols, types):
        return pa.table({c: pa.array([r[i] for r in rows], type=t) for i, (c, t) in enumerate(zip(cols, types))})

    con.register("_snp", _arrow(snp_rows, ["rsid", "gene", "chrom", "pos", "orientation", "summary", "gmaf", "text"],
                                [pa.string(), pa.string(), pa.string(), pa.int64(), pa.string(), pa.string(), pa.float64(), pa.string()]))
    con.execute("DROP TABLE IF EXISTS snpedia_snp")
    con.execute("CREATE TABLE snpedia_snp AS SELECT * FROM _snp")
    con.register("_geno", _arrow(geno_rows, ["rsid", "allele1", "allele2", "magnitude", "repute", "summary"],
                                 [pa.string(), pa.string(), pa.string(), pa.float64(), pa.string(), pa.string()]))
    con.execute("DROP TABLE IF EXISTS snpedia_geno")
    con.execute("CREATE TABLE snpedia_geno AS SELECT * FROM _geno")
    con.register("_alias", _arrow(alias_rows, ["iid", "rsid"], [pa.string(), pa.string()]))
    con.execute("DROP TABLE IF EXISTS snpedia_alias")
    con.execute("CREATE TABLE snpedia_alias AS SELECT * FROM _alias")
    for v in ("_snp", "_geno", "_alias"):
        con.unregister(v)
    con.execute("CREATE INDEX IF NOT EXISTS snpedia_geno_rsid ON snpedia_geno(rsid)")
    con.execute("CREATE INDEX IF NOT EXISTS snpedia_snp_rsid ON snpedia_snp(rsid)")
    console.print(f"  snpedia_snp: {len(snp_rows):,} · snpedia_geno: {len(geno_rows):,} · alias 23andMe i-ids: {len(alias_rows):,}")


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ----------------------------------------------------------------------------- ClinVar / GWAS / DGIdb
def build_clinvar(con: duckdb.DuckDBPyConnection, path: Path) -> None:
    con.execute("DROP TABLE IF EXISTS clinvar")
    con.execute(f"""
        CREATE TABLE clinvar AS
        SELECT 'rs' || "RS# (dbSNP)" AS rsid,
               GeneSymbol AS gene, Name AS name, Type AS vtype,
               ClinicalSignificance AS clinsig, ReviewStatus AS review_status,
               CASE ReviewStatus
                    WHEN 'practice guideline' THEN 4
                    WHEN 'reviewed by expert panel' THEN 3
                    WHEN 'criteria provided, multiple submitters, no conflicts' THEN 2
                    WHEN 'criteria provided, single submitter' THEN 1
                    WHEN 'criteria provided, conflicting classifications' THEN 1
                    WHEN 'criteria provided, conflicting interpretations' THEN 1
                    ELSE 0 END AS stars,
               PhenotypeList AS phenotypes, Chromosome AS chrom, TRY_CAST(Start AS INTEGER) AS pos,
               ReferenceAlleleVCF AS ref, AlternateAlleleVCF AS alt, VariationID AS variation_id,
               TRY_CAST(NumberSubmitters AS INTEGER) AS n_submitters, LastEvaluated AS last_evaluated
        FROM read_csv('{path.as_posix()}', delim='\t', header=true, all_varchar=true, quote='', escape='', ignore_errors=true)
        WHERE Assembly = 'GRCh37' AND "RS# (dbSNP)" <> '-1'
          AND (ClinicalSignificance ILIKE '%pathogenic%' OR ClinicalSignificance ILIKE '%risk factor%'
               OR ClinicalSignificance ILIKE '%drug response%' OR ClinicalSignificance ILIKE '%protective%'
               OR ClinicalSignificance ILIKE '%affects%')
    """)
    con.execute("CREATE INDEX IF NOT EXISTS clinvar_rsid ON clinvar(rsid)")
    console.print(f"  clinvar: {con.execute('select count(*) from clinvar').fetchone()[0]:,} filas")


def build_gwas(con: duckdb.DuckDBPyConnection, zpath: Path) -> None:
    with zipfile.ZipFile(zpath) as z:
        name = z.namelist()[0]
        tsv = RAW_DIR / "gwas_associations_full.tsv"
        if not tsv.exists():
            with z.open(name) as src, open(tsv, "wb") as dst:
                while chunk := src.read(1 << 24):
                    dst.write(chunk)
    con.execute("DROP TABLE IF EXISTS gwas")
    con.execute(f"""
        CREATE TABLE gwas AS
        SELECT lower(trim(split_part("STRONGEST SNP-RISK ALLELE", '-', 1))) AS rsid,
               upper(trim(split_part("STRONGEST SNP-RISK ALLELE", '-', 2))) AS risk_allele,
               "DISEASE/TRAIT" AS trait, MAPPED_TRAIT AS mapped_trait,
               TRY_CAST("P-VALUE" AS DOUBLE) AS pvalue,
               TRY_CAST("OR or BETA" AS DOUBLE) AS or_beta,
               "95% CI (TEXT)" AS ci_text,
               CASE WHEN "95% CI (TEXT)" ILIKE '%increase%' OR "95% CI (TEXT)" ILIKE '%decrease%' OR "95% CI (TEXT)" ILIKE '%unit%' OR "95% CI (TEXT)" ILIKE '%SD%'
                    THEN 'beta' ELSE 'or' END AS effect_type,
               COALESCE(NULLIF(MAPPED_GENE, ''), "REPORTED GENE(S)") AS gene,
               PUBMEDID AS pubmedid, "FIRST AUTHOR" AS first_author, "DATE" AS pub_date,
               "INITIAL SAMPLE SIZE" AS sample_text,
               TRY_CAST("RISK ALLELE FREQUENCY" AS DOUBLE) AS raf,
               "STUDY ACCESSION" AS study
        FROM read_csv('{tsv.as_posix()}', delim='\t', header=true, all_varchar=true, quote='', escape='', ignore_errors=true)
        WHERE "STRONGEST SNP-RISK ALLELE" LIKE 'rs%' AND "STRONGEST SNP-RISK ALLELE" NOT LIKE '%;%'
          AND TRY_CAST("P-VALUE" AS DOUBLE) <= 5e-8
    """)
    con.execute("CREATE INDEX IF NOT EXISTS gwas_rsid ON gwas(rsid)")
    tsv.unlink(missing_ok=True)  # keep only the zip in raw/ to respect the budget
    console.print(f"  gwas: {con.execute('select count(*) from gwas').fetchone()[0]:,} asociaciones (p<=5e-8)")


def build_dgidb(con: duckdb.DuckDBPyConnection, path: Path) -> None:
    con.execute("DROP TABLE IF EXISTS dgidb")
    con.execute(f"""
        CREATE TABLE dgidb AS
        SELECT gene_name AS gene, drug_name AS drug, NULLIF(interaction_type,'NULL') AS interaction_type,
               TRY_CAST(interaction_score AS DOUBLE) AS score, approved = 'TRUE' AS approved,
               interaction_source_db_name AS source
        FROM read_csv('{path.as_posix()}', delim='\t', header=true, all_varchar=true, quote='', escape='', ignore_errors=true)
        WHERE gene_name IS NOT NULL AND drug_name IS NOT NULL
    """)
    console.print(f"  dgidb: {con.execute('select count(*) from dgidb').fetchone()[0]:,} interacciones")


# ----------------------------------------------------------------------------- Open Targets
def fetch_opentargets(client: httpx.Client) -> None:
    for ds in OT_DATASETS:
        d = RAW_DIR / "opentargets" / ds
        d.mkdir(parents=True, exist_ok=True)
        if any(d.glob("*.parquet")):
            continue
        listing = client.get(f"{OT_BASE}/{ds}/", timeout=120).text
        files = re.findall(r'href="([^"]+\.parquet)"', listing)
        if not files:
            raise SystemExit(f"Open Targets: no se encontraron parquet en {ds}")
        for f in files:
            download(f"{OT_BASE}/{ds}/{f}", d / f, client, label=f"OT {ds}")


def build_opentargets(con: duckdb.DuckDBPyConnection) -> None:
    base = (RAW_DIR / "opentargets").as_posix()
    con.execute("DROP TABLE IF EXISTS ot_molecule")
    con.execute(f"""
        CREATE TABLE ot_molecule AS
        SELECT id AS chembl_id, name, drugType AS drug_type, maximumClinicalStage AS max_stage, description,
               list_transform(synonyms, x -> x.label) AS synonyms, list_transform(tradeNames, x -> x.label) AS trade_names
        FROM read_parquet('{base}/drug_molecule/*.parquet')""")
    con.execute("DROP TABLE IF EXISTS ot_moa")
    con.execute(f"""
        CREATE TABLE ot_moa AS
        SELECT unnest(chemblIds) AS chembl_id, actionType AS action_type, mechanismOfAction AS mechanism,
               targetName AS target_name, targets AS target_ids
        FROM read_parquet('{base}/drug_mechanism_of_action/*.parquet')""")
    con.execute("DROP TABLE IF EXISTS ot_target")
    con.execute(f"""
        CREATE TABLE ot_target AS
        SELECT id AS ensg, approvedSymbol AS symbol, approvedName AS name, biotype,
               array_to_string(functionDescriptions, ' ') AS function_text,
               list_transform(pathways, x -> x.pathway) AS pathways,
               list_transform(list_filter(tractability, x -> x.value), x -> x.modality || ':' || x.id) AS tractability
        FROM read_parquet('{base}/target/*.parquet')
        WHERE biotype = 'protein_coding' OR approvedSymbol LIKE 'MT-%' OR approvedSymbol LIKE 'MIR%' """)
    con.execute("DROP TABLE IF EXISTS ot_clinical_target")
    con.execute(f"""
        CREATE TABLE ot_clinical_target AS
        SELECT drugId AS chembl_id, targetId AS ensg, maxClinicalStage AS max_stage,
               list_transform(diseases, x -> x.diseaseFromSource)[1:8] AS diseases
        FROM read_parquet('{base}/clinical_target/*.parquet')""")
    con.execute("DROP TABLE IF EXISTS ot_pgx")
    con.execute(f"""
        CREATE TABLE ot_pgx AS
        SELECT lower(variantRsId) AS rsid, genotype, genotypeAnnotationText AS annotation, evidenceLevel AS level,
               pgxCategory AS category, phenotypeText AS phenotype, targetFromSourceId AS ensg,
               list_transform(drugs, x -> x.drugFromSource) AS drugs, literature, datasourceId AS source,
               isDirectTarget AS direct_target,
               list_transform(variantAnnotation, x -> x.comparisonAlleleOrGenotype) AS comparisons,
               list_transform(variantAnnotation, x -> x.baseAlleleOrGenotype) AS bases
        FROM read_parquet('{base}/pharmacogenomics/*.parquet')
        WHERE variantRsId IS NOT NULL""")
    con.execute("CREATE INDEX IF NOT EXISTS ot_pgx_rsid ON ot_pgx(rsid)")
    for t in ["ot_molecule", "ot_moa", "ot_target", "ot_clinical_target", "ot_pgx"]:
        console.print(f"  {t}: {con.execute(f'select count(*) from {t}').fetchone()[0]:,}")


# ----------------------------------------------------------------------------- CPIC
def fetch_cpic(client: httpx.Client) -> None:
    d = RAW_DIR / "cpic"
    d.mkdir(parents=True, exist_ok=True)
    endpoints = {
        "pair.json": "pair?select=genesymbol,drugid,cpiclevel,clinpgxlevel,pgxtesting,citations,guidelineid&removed=eq.false",
        "drug.json": "drug?select=drugid,name,atcid",
        "allele.json": "allele?select=genesymbol,name,clinicalfunctionalstatus,activityvalue,strength,findings",
        "gene.json": "gene?select=symbol,chr,lookupmethod",
        "guideline.json": "guideline?select=id,name,url,genes",
    }
    for fname, q in endpoints.items():
        p = d / fname
        if p.exists():
            continue
        r = client.get(f"{CPIC_API}/{q}", headers={"Prefer": "count=none"}, timeout=120)
        r.raise_for_status()
        p.write_text(r.text, encoding="utf-8")


def build_cpic(con: duckdb.DuckDBPyConnection) -> None:
    d = (RAW_DIR / "cpic").as_posix()
    con.execute("DROP TABLE IF EXISTS cpic_pair")
    con.execute(f"""
        CREATE TABLE cpic_pair AS
        SELECT p.genesymbol AS gene, dr.name AS drug, p.cpiclevel AS cpic_level, p.clinpgxlevel AS pgkb_level,
               p.pgxtesting AS pgx_testing, p.citations, g.url AS guideline_url
        FROM read_json_auto('{d}/pair.json') p
        LEFT JOIN read_json_auto('{d}/drug.json') dr ON dr.drugid = p.drugid
        LEFT JOIN read_json_auto('{d}/guideline.json') g ON g.id = p.guidelineid""")
    con.execute("DROP TABLE IF EXISTS cpic_allele")
    con.execute(f"""
        CREATE TABLE cpic_allele AS
        SELECT genesymbol AS gene, name AS allele, clinicalfunctionalstatus AS function, activityvalue AS activity,
               strength, findings
        FROM read_json_auto('{d}/allele.json')""")
    console.print(f"  cpic_pair: {con.execute('select count(*) from cpic_pair').fetchone()[0]:,} · cpic_allele: {con.execute('select count(*) from cpic_allele').fetchone()[0]:,}")


# ----------------------------------------------------------------------------- main
def fetch_all(skip_snpedia: bool = False) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    m = _manifest()
    with _http() as client:
        console.rule("[bold]1/3 Descargas")
        for key, (url, fname) in FILES.items():
            p = download(url, RAW_DIR / fname, client, label=key)
            m["sources"][key] = {"url": url, "bytes": p.stat().st_size, "fetched": _now()}
            _budget_check()
        fetch_opentargets(client)
        m["sources"]["opentargets"] = {"url": OT_BASE, "datasets": OT_DATASETS, "fetched": _now()}
        fetch_cpic(client)
        m["sources"]["cpic"] = {"url": CPIC_API, "fetched": _now()}
        _budget_check()
        if not skip_snpedia:
            fetch_snpedia(client)
            m["sources"]["snpedia"] = {"url": SNPEDIA_API, "fetched": _now(), "license": "CC-BY-NC-SA 3.0"}
        _budget_check()
    _save_manifest(m)

    console.rule("[bold]2/3 Construcción de la base DuckDB")
    tmp = DB_PATH.with_suffix(".building.duckdb")
    tmp.unlink(missing_ok=True)
    con = duckdb.connect(str(tmp))
    try:
        build_clinvar(con, RAW_DIR / FILES["clinvar"][1])
        build_gwas(con, RAW_DIR / FILES["gwas"][1])
        build_dgidb(con, RAW_DIR / FILES["dgidb"][1])
        build_opentargets(con)
        build_cpic(con)
        if (RAW_DIR / "snpedia" / "snps.jsonl").exists():
            build_snpedia_tables(con)
        con.execute("CREATE TABLE IF NOT EXISTS meta(key VARCHAR, value VARCHAR)")
        con.execute("INSERT INTO meta VALUES ('built', ?)", [_now()])
        con.execute("CHECKPOINT")
    finally:
        con.close()
    if DB_PATH.exists():
        DB_PATH.unlink()
    tmp.rename(DB_PATH)

    console.rule("[bold]3/3 Resumen")
    raw = _raw_bytes()
    console.print(f"Descargado: [bold]{raw/1e9:.2f} GB[/] (presupuesto 2 GB) · DB: [bold]{DB_PATH.stat().st_size/1e6:.0f} MB[/] → {DB_PATH}")
    m["built"] = _now()
    m["raw_bytes"] = raw
    m["db_bytes"] = DB_PATH.stat().st_size
    _save_manifest(m)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")
