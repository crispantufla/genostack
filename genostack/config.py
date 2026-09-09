"""Paths, constants and significance thresholds."""
from __future__ import annotations

import os
from pathlib import Path

PKG_DIR = Path(__file__).resolve().parent
PANEL_DIR = PKG_DIR / "panel"

DATA_DIR = Path(os.environ.get("GENOSTACK_DATA", Path.cwd() / "data")).resolve()
RAW_DIR = DATA_DIR / "raw"
DB_PATH = DATA_DIR / "genostack.duckdb"
MANIFEST_PATH = DATA_DIR / "manifest.json"

DOWNLOAD_BUDGET_BYTES = 2 * 1024**3  # hard budget for pre-downloaded data

# ---- significance thresholds (see PLAN.md §4) ----
SNPEDIA_MIN_MAGNITUDE = 2.5
SNPEDIA_HIGHLIGHT_MAGNITUDE = 3.0
CLINVAR_MIN_STARS = 2
GWAS_MAX_P = 5e-8
GWAS_MIN_OR = 1.3
GWAS_MAX_PLAUSIBLE_OR = 10.0      # above this it is a catalog artefact, not biology
GWAS_STRONG_OR = 5.0              # 5-10 is possible but demands a well-powered study
GWAS_STRONG_OR_MIN_N = 20_000
GWAS_MIN_BETA_SD = 0.10
PANEL_MIN_IMPACT = 1

# Evidence-level weights used by the composite score
EVIDENCE_WEIGHT = {"A": 1.0, "B": 0.8, "C": 0.55, "D": 0.35, "E": 0.2}

# GWAS trait keyword whitelist (lower-case substrings) — categories relevant to the goals
GWAS_TRAIT_KEYWORDS = [
    # nutrients / biomarkers
    "vitamin", "folate", "b12", "cobalamin", "homocysteine", "iron", "ferritin", "transferrin",
    "hemoglobin", "haemoglobin", "zinc", "magnesium", "selenium", "copper", "calcium", "phosphate",
    "urate", "uric acid", "omega", "docosahexaenoic", "eicosapentaenoic", "fatty acid", "choline",
    "carnitine", "coenzyme", "glutathione", "25-hydroxy", "retinol", "carotene", "tocopherol",
    # metabolic / cardio
    "glucose", "insulin", "hba1c", "glycated", "type 2 diabetes", "triglyceride", "cholesterol",
    "ldl", "hdl", "apolipoprotein", "lipoprotein", "blood pressure", "hypertension",
    "coronary", "myocardial", "atrial fibrillation", "stroke", "thrombo", "venous",
    "body mass", "obesity", "waist", "adiposity", "fat mass", "lean mass", "liver enzyme",
    "alanine aminotransferase", "gamma-glutamyl", "fatty liver", "nafld", "gout",
    # cognition / neuro / mood
    "cognitive", "intelligence", "memory", "educational attainment", "reaction time",
    "alzheimer", "dementia", "parkinson", "depress", "anxiety", "neurotic", "adhd", "attention",
    "bipolar", "schizophrenia", "autism", "migraine", "sleep", "insomnia", "chronotype",
    "caffeine", "coffee", "alcohol", "smoking", "nicotine",
    # inflammation / immune
    "c-reactive", "crp", "interleukin", "inflammat", "autoimmun", "celiac", "coeliac", "crohn",
    "colitis", "rheumatoid", "psoriasis", "asthma", "allerg", "eczema", "lupus", "thyroid",
    "hypothyroid", "tsh", "vitiligo", "multiple sclerosis",
    # longevity / performance / misc
    "longevity", "lifespan", "parental", "telomere", "testosterone", "estradiol", "shbg",
    "igf", "growth hormone", "bone mineral", "osteopor", "fracture", "muscle", "grip strength",
    "vo2", "athlet", "endurance", "hair loss", "baldness", "macular", "glaucoma", "hearing",
    "kidney", "egfr", "creatinine", "cystatin", "albumin", "prostate", "breast cancer",
    "colorectal cancer", "melanoma", "lung cancer", "cancer",
]

GOALS = ["health", "energy", "cognition", "longevity", "performance", "pharmaco"]
