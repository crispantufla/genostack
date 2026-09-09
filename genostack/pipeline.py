"""End-to-end `analyze` pipeline: parse → annotate (local DB + panel) → rank → research (ACP agent) → report."""
from __future__ import annotations

import asyncio
import time
from pathlib import Path

from rich.console import Console
from rich.live import Live
from rich.table import Table

from .acp_client import AgentRunner
from .annotate import annotate, gene_context, open_db
from .evidence import ensembl_frequencies
from .panel import load_panel
from .report import print_report, render_markdown
from .research import run_research
from .score import rank

console = Console()


def _add_alphagenome(analysis, console) -> None:
    """Predice el efecto molecular de los hallazgos cuyo mecanismo no se conoce ya.

    Se anota en `extra` y viaja al agente y al informe SIEMPRE etiquetado como predicción: no altera
    el score ni el grado de evidencia de nada.
    """
    from .alphagenome import predict
    cand = [(f.rsid, f.genotype) for b in analysis.modules for f in b.findings
            if f.source in ("gwas", "snpedia") and "," not in f.rsid and len(f.genotype) == 2]
    cand = list(dict.fromkeys(cand))
    if not cand:
        console.print("  AlphaGenome: ningún hallazgo candidato (todos con mecanismo conocido)")
        return
    console.rule("[bold]Predicción molecular (AlphaGenome)")
    console.print(f"  {len(cand)} variantes candidatas · predicción computacional, no observación")
    try:
        preds = predict(cand, progress=lambda n, tot, p: console.print(
            f"    [{n}/{tot}] {p.rsid}: {p.summary(2)}"))
    except RuntimeError as e:
        console.print(f"[yellow]  AlphaGenome no disponible: {e!s:.180}[/]")
        return
    hit = 0
    for b in analysis.modules:
        for f in b.findings:
            p = preds.get(f.rsid)
            if p and p.effects:
                f.extra["alphagenome"] = {"summary": p.summary(), "effects": p.effects[:4],
                                          "consequence": p.consequence}
                hit += 1
    console.print(f"  {hit} hallazgos con predicción molecular añadida")


def _profile_summary(analysis) -> str:
    """Compact whole-genome summary shared with every module prompt (so modules can see cross-module context)."""
    lines = []
    for b in analysis.modules:
        if b.id == "other":
            continue
        items = [f"{f.gene} {f.genotype} ({f.detail[:60].strip()})" for f in b.findings[:8]]
        if items:
            lines.append(f"- {b.title}: " + "; ".join(items))
    other = [b for b in analysis.modules if b.id == "other"]
    if other and other[0].findings:
        lines.append("- Otros (auto): " + "; ".join(f"{f.title[:50]} {f.genotype}" for f in other[0].findings[:10]))
    return "\n".join(lines)


async def run_analysis(file: Path, agent_cmd: str | None, concurrency: int = 3, lang: str = "es", save: bool = True,
                       goals: list[str] | None = None, use_alphagenome: bool = False) -> None:
    goals = goals or ["health", "energy", "cognition", "longevity"]
    t0 = time.time()
    from .parsers import parse_raw
    console.rule("[bold]GENOSTACK")
    console.print(f"Leyendo [bold]{file}[/]…")
    genome = parse_raw(file)
    console.print(f"  {genome.source_format} · {genome.build} · {genome.n_calls:,} genotipos ({genome.n_nocalls:,} no-calls) · sexo inferido: {genome.sex_guess()}")
    if genome.build != "GRCh37":
        console.print("[yellow]  Aviso: el archivo no parece GRCh37; los rsIDs siguen siendo válidos, pero las posiciones no se usan.[/]")

    con = open_db()
    try:
        db_built = con.execute("SELECT value FROM meta WHERE key='built'").fetchone()
        console.print("Cruzando con bases locales (SNPedia, ClinVar, GWAS, PharmGKB) y panel curado…")
        findings, stats = annotate(genome, con)
        panel = load_panel()
        analysis = rank(findings, panel, goals)
        analysis.stats = stats
        n_major = sum(len(b.findings) for b in analysis.modules)
        console.print(f"  {len(findings):,} hallazgos crudos → {n_major} significativos en {len(analysis.modules)} módulos")

        # gene context per module for the agent
        gene_ctx_by_module = {}
        if agent_cmd:
            for b in analysis.modules:
                if b.id != "other":
                    gene_ctx_by_module[b.id] = gene_context(con, b.genes[:20])
    finally:
        con.close()

    if use_alphagenome:
        _add_alphagenome(analysis, console)

    # online: allele frequencies for the prioritized rsids (cheap, informative)
    research = synth = None
    agent_log_dir = None
    agent_name = None
    if agent_cmd:
        rsids = [f.rsid for b in analysis.modules for f in b.findings if f.source != "haplotype"][:200]
        try:
            freqs = await asyncio.wait_for(ensembl_frequencies(rsids), timeout=120)
            for b in analysis.modules:
                for f in b.all_findings():
                    if f.rsid in freqs:
                        f.extra["freq"] = freqs[f.rsid]
        except Exception:  # noqa: BLE001
            pass

        console.rule("[bold]Investigación con agente ACP")
        status: dict[str, str] = {}

        def make_table() -> Table:
            t = Table(show_header=True, header_style="bold", expand=False)
            t.add_column("Tarea")
            t.add_column("Estado")
            for k, v in status.items():
                t.add_row(k, v)
            return t

        def progress(label: str, text: str) -> None:
            status[label] = text

        profile = _profile_summary(analysis)
        try:
            async with AgentRunner(agent_cmd, progress=progress, concurrency=concurrency) as runner:
                agent_log_dir = runner.log_dir
                agent_name = runner.agent_info
                console.print(f"  Agente: {runner.agent_info} · {len([b for b in analysis.modules if b.id != 'other'])} módulos · concurrencia {concurrency}")
                with Live(make_table(), console=console, refresh_per_second=2, transient=True) as live:
                    async def ticker():
                        while True:
                            live.update(make_table())
                            await asyncio.sleep(0.5)
                    tk = asyncio.create_task(ticker())
                    try:
                        research, audits, synth = await run_research(runner, analysis, gene_ctx_by_module, profile, lang, goals,
                                                                     log=lambda m: console.print("  " + m))
                    finally:
                        tk.cancel()
                ok = sum(1 for r in research if not r.error)
                console.print(f"  Módulos investigados: {ok}/{len(research)} · síntesis: {'ok' if synth and not synth.error else 'fallida'}")
        except Exception as e:  # noqa: BLE001
            console.print(f"[red]No se pudo usar el agente ACP ({type(e).__name__}: {e!s:.200}). Se imprime el análisis local.[/]")
            research = research or None

    meta = {"file": str(file), "format": genome.source_format, "build": genome.build, "sex": genome.sex_guess(),
            "agent": agent_name, "db_built": db_built[0] if db_built else "?", "agent_log_dir": str(agent_log_dir) if agent_log_dir else None}
    md = render_markdown(analysis, research, synth, stats, meta, lang)
    console.rule("[bold]Informe")
    print_report(md, console)
    if save:
        out = file.with_suffix(file.suffix + ".genostack.md")
        out.write_text(md, encoding="utf-8")
        console.print(f"\n[green]Informe guardado en[/] {out}")
    console.print(f"[dim]Tiempo total: {time.time()-t0:.0f}s[/]")
