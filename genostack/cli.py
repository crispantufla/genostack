"""genostack CLI: `fetch-data` and `analyze`."""
from __future__ import annotations

import argparse
import asyncio
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> None:
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8")
        except Exception:
            pass
    ap = argparse.ArgumentParser(prog="genostack", description="Análisis genético local + investigación de intervenciones vía agente ACP")
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch-data", help="Descarga e indexa las bases locales (≤ 2 GB). Reanudable.")
    f.add_argument("--skip-snpedia", action="store_true", help="No descargar SNPedia (más rápido; menos cobertura)")

    a = sub.add_parser("analyze", help="Analiza un raw report y escribe el informe en consola")
    a.add_argument("file", help="Raw report (23andMe/Ancestry/MyHeritage/FTDNA/VCF; .txt/.csv/.zip/.gz)")
    a.add_argument("--agent", default="npx -y @agentclientprotocol/claude-agent-acp",
                   help="Comando del agente ACP (por defecto Claude Code vía claude-agent-acp)")
    a.add_argument("--offline", action="store_true", help="No usar agente ni red: sólo análisis local")
    a.add_argument("--concurrency", type=int, default=3, help="Sesiones ACP simultáneas")
    a.add_argument("--lang", default="es", choices=["es", "en"], help="Idioma del informe")
    a.add_argument("--no-save", action="store_true", help="No guardar copia .md junto al input")
    a.add_argument("--goals", default="health,energy,cognition,longevity",
                   help="Objetivos a ponderar (health,energy,cognition,longevity,performance,pharmaco)")

    args = ap.parse_args(argv)
    if args.cmd == "fetch-data":
        from .fetch import fetch_all
        fetch_all(skip_snpedia=args.skip_snpedia)
    elif args.cmd == "analyze":
        from .pipeline import run_analysis
        asyncio.run(run_analysis(Path(args.file), agent_cmd=None if args.offline else args.agent,
                                 concurrency=args.concurrency, lang=args.lang, save=not args.no_save,
                                 goals=[g.strip() for g in args.goals.split(",") if g.strip()]))


if __name__ == "__main__":
    main()
