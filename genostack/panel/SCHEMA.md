# Panel curado — esquema YAML

Un fichero por módulo funcional: `genostack/panel/<module>.yaml`.

**Regla de oro de alelos:** TODOS los alelos se escriben en la orientación de la hebra + de GRCh37, tal y como aparecen en un raw file de 23andMe/Ancestry (es decir, como figuran en dbSNP "forward"/ SNPedia con orientation=plus). Si SNPedia lista el SNP con `orientation=minus`, los alelos de SNPedia deben complementarse (A↔T, C↔G) antes de escribirlos aquí. El genotipo del usuario se compara contra `genotypes` tras ordenar alfabéticamente sus dos letras (p.ej. "GA" → "AG"). Para cromosoma X en varones / MT / Y puede haber un único alelo: incluir también claves de una letra cuando sea relevante.

```yaml
module: methylation                 # id corto, snake_case, único
title: "Metilación / ciclo de un carbono"
goals: [health, energy, cognition]  # de: health, energy, cognition, longevity, performance, pharmaco
description: >
  Una o dos frases sobre qué cubre el módulo y por qué importa.

variants:
  - id: MTHFR_C677T                  # único en todo el panel
    gene: MTHFR
    rsid: rs1801133
    name: "C677T (Ala222Val)"
    effect_allele: A                 # alelo de efecto, hebra + GRCh37
    other_allele: G
    direction: risk                  # risk | beneficial | mixed  (efecto del alelo de efecto)
    genotypes:                       # claves: 2 letras ordenadas alfabéticamente; impact 0-4
      AA: {phenotype: "Homocigoto 677TT: actividad MTHFR ~30 %, homocisteína ↑ si folato/B2 bajos", impact: 3}
      AG: {phenotype: "Heterocigoto: actividad ~65 %", impact: 1}
      GG: {phenotype: "Normal", impact: 0}
    effect_size: "Actividad enzimática 30 % (TT) / 65 % (CT) vs 100 %; TT: homocisteína +~25 %"
    evidence: A                      # A meta-análisis/RCT/consenso; B replicado; C limitado; D mecanístico/preliminar
    mechanism: "Variante termolábil; menor 5-MTHF → menor remetilación de homocisteína y menos SAMe"
    levers:                          # palancas semilla (el agente las verificará y ampliará)
      - "L-metilfolato 400-1000 µg/d en lugar de ácido fólico"
      - "Riboflavina 1.6 mg/d: reduce homocisteína ~40 % en TT (RCT)"
      - "Medir homocisteína, folato eritrocitario, B12"
    caveats: "Efecto clínico depende del estado de folato; en homocigotos, evitar altas dosis de ácido fólico sintético"
    refs: ["PMID:9545395", "PMID:20101055"]
```

Campos obligatorios: `id, gene, rsid, effect_allele, other_allele, direction, genotypes, effect_size, evidence, mechanism, levers`. `impact` obligatorio en cada genotipo (0 = sin efecto; 1 leve; 2 moderado; 3 fuerte; 4 clínico/muy fuerte). Sólo genotipos con impact ≥ 1 acaban en el informe, así que sé honesto: efectos marginales → impact 0 o 1, no 2-3.

## Haplotipos / combinaciones (`haplotypes.yaml`)

```yaml
rules:
  - id: APOE
    gene: APOE
    title: "Genotipo APOE"
    type: combo                      # combo = coincidencia exacta de genotipos
    rsids: [rs429358, rs7412]
    combos:
      "rs429358=TT;rs7412=CC": {name: "ε3/ε3", phenotype: "Referencia", impact: 0}
      "rs429358=CT;rs7412=CC": {name: "ε3/ε4", phenotype: "Riesgo Alzheimer ×3, LDL ↑", impact: 3, direction: risk}
      "rs429358=CC;rs7412=CC": {name: "ε4/ε4", phenotype: "Riesgo Alzheimer ×12-15", impact: 4, direction: risk}
      "rs429358=TT;rs7412=CT": {name: "ε2/ε3", phenotype: "LDL ↓, Alzheimer ↓; riesgo disbetalipoproteinemia si ε2/ε2", impact: 1, direction: beneficial}
      "rs429358=TT;rs7412=TT": {name: "ε2/ε2", phenotype: "...", impact: 2, direction: mixed}
      "rs429358=CT;rs7412=CT": {name: "ε2/ε4 (o ε1/ε3, raro)", phenotype: "...", impact: 2, direction: mixed}
    evidence: A
    mechanism: "..."
    levers: ["..."]
    refs: ["PMID:..."]

  - id: CYP2C19
    gene: CYP2C19
    title: "Fenotipo metabolizador CYP2C19"
    type: allele_count               # suma de alelos de pérdida/ganancia de función
    alleles:
      - {star: "*2",  rsid: rs4244285,  allele: A, function: loss}
      - {star: "*3",  rsid: rs4986893,  allele: A, function: loss}
      - {star: "*17", rsid: rs12248560, allele: T, function: gain}
    phenotypes:                      # evaluadas en orden; primera que cumple gana
      - {when: "loss>=2", name: "Metabolizador lento (PM)", impact: 4, direction: risk}
      - {when: "loss==1", name: "Metabolizador intermedio (IM)", impact: 2, direction: risk}
      - {when: "gain>=2", name: "Metabolizador ultrarrápido (UM)", impact: 3, direction: mixed}
      - {when: "gain==1", name: "Metabolizador rápido (RM)", impact: 2, direction: mixed}
      - {when: "true",   name: "Metabolizador normal (NM)", impact: 0}
    evidence: A
    mechanism: "..."
    levers: ["..."]
    refs: ["PMID:..."]
```

Las condiciones `when` sólo usan `loss`, `gain`, `count` (nº de alelos de efecto en total), comparadores `==, >=, <=, >, <` y `true`. Si falta algún rsid en el archivo del usuario, la regla se evalúa con los disponibles y se marca `partial: true` en el resultado.
