# Pangenome-Scale Reconstruction of Lactobacillaceae Metabolism
Automated reconstruction of strain-specific Genome-Scale Metabolic Models (GEMs) across the Lactobacillaceae family, with gap filling and downstream comparative analyses.

> End-to-end workflow: curated inputs → GEM generation → gap filling → QC/curation → flux/essentiality analyses → niche enrichment and panGEM FVA reports.

---

## Table of Contents
- [Overview](#overview)
- [Workflow](#workflow)
  - [1) GEM Generation (`GEMgenerator.py`)](#1-gem-generation-gemgeneratorpy)
  - [2) Gap Filling (`dgap.py`, `run_dgap.py`)](#2-gap-filling-dgappy-rundgappy)
  - [3) Manual Curation](#3-manual-curation)
  - [4) Reporting & Analysis Tools](#4-reporting--analysis-tools)
- [Figure](#figure)
- [Directory Layout](#directory-layout)
- [Dependencies](#dependencies)
- [Quick Start](#quick-start)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [License](#license)

---

## Overview
This repository provides code and documentation to reconstruct **strain-specific GEMs** for Lactobacillaceae. The workflow is optimized for high-performance environments (e.g., 96 cores on Azure) but can be adapted to smaller or larger machines.


<p align="center">
  <img src="docs/msystems.00156-24.f001.jpg" width="450" alt="Automated pipeline: genomes → QA/QC → annotation (BAKTA) → pangenome (CD-HIT) → GEM reconstruction (CarveMe) → panGPRs → neighborhood analysis; with protein stoichiometry, 3D modeling, and structural analysis integration.">
</p>

<p align="center"><em>Figure 1. High-level workflow for automated pangenome and metabolic modeling.</em></p>


---

## Workflow

### 1) GEM Generation (`GEMgenerator.py`)
Generates draft strain-specific GEMs.

**Inputs**
- Reactome model (JSON)
- Gene sequences for Reactome GPRs:
  - Nucleotide FASTA (`.fna`)
  - Amino acid FASTA (`.faa`)
- GenBank files (`.gbk`) for target strains

**Output**
- Draft strain-specific GEMs (one model per strain)

**Example**
```bash
python GEMgenerator.py \
  --reactome data/reactome_model.json \
  --genes-fna data/sequences/genes.fna \
  --genes-faa data/sequences/genes.faa \
  --genbank  data/genomes/ \
  --out      results/gems_draft \
  --threads  96
```

## Additional summary (03/24)
### Core Pipeline

GEMgenerator.py builds draft strain-specific GEMs from a Reactome reference model plus genome files. It parses GenBank annotations, runs BLAST-based orthology searches, constructs presence/absence matrices, prunes the reference model per strain, renames genes to strain-specific loci, and writes draft JSON models.

dgap.py is the main utility library for gapfilling and most downstream model analysis. It defines the m9() media setup, scans model directories for feasible vs failed models, identifies candidate gapfilled reactions from a template model, tests whether reactions are essential, and also includes many later helper functions for exchange fluxes, essential reactions, carbon-source tests, product-associated reaction knockouts, connectivity, and model cleanup. It is the real “toolbox” file in this repo.

run_dgap.py is a hard-coded driver script for one particular gapfilling run. It calls scan, tempfind, gaps, zx, and fluxanalyze from dgap.py, then adds selected reactions from one template model into failed models and splits outputs into feasible/ vs failed2/.

dgap2.py is a variant of dgap.py. It overlaps heavily, but uses different media bounds and a few extra utility functions near the bottom. It looks like an experimental fork rather than a clean successor.

### Downstream Analysis Scripts

basic_inf.py loops over genus folders under /home/omidard/allgems, reads corrected models from each corr/ directory, applies m9(), and writes a basic_inf.csv with model id, growth rate, number of reactions, number of gap reactions, and number of genes.

basic_counts.py parallelizes dgap.basic_counts() over /home/omidard/gems/allgems/ and writes one summary counts table.

species_all_gaps.py collects all reactions associated with the artificial GAP gene from each corrected model and writes a genus-level presence matrix of gapfilled reactions.

species_all_fluxes.py loads every model in /home/omidard/gems/allgems/, applies m9(), runs FBA, and writes a matrix of all reaction fluxes across all species.

species_knockout_fluxes.py tries to assess how reaction knockouts affect exchange fluxes and then plot the distributions. It operates per genus folder and writes species_knockout_exchange_fluxes.csv plus a figure. This one looks especially brittle because of how it reloads and mutates models inside nested loops.

growth_rates.py augments an existing metadata CSV with wild-type biomass flux and a list of exchange reaction fluxes for each model.

coreflux.py parallelizes dgap.coreflux() to build a cross-model table of only the nonzero reaction fluxes.

lacto_pangem_fva.py is the cleanest standalone analysis script in the repo. It runs full FVA and parsimonious FBA in parallel across JSON models in /home/omidard/PanGEMs, then writes lactopan_fva.csv and lactopan_parsFBA.csv.

### Carbon / Product / Niche Scripts

carb_sources.py tests growth on a predefined list of carbon sources by turning glucose off, enabling one carbon exchange at a time, and recording BIOMASS2 flux per model.

carb.py is a post-processing script. It reads earlier CSV outputs like minimal_media.csv and carbon_sources.csv, reshapes them, and combines them into one larger table inf.csv.

pr_ess.py finds producer models for a target product, computes active reactions, screens nonessential knockouts, and writes a product-based single-knockout table. In the current code it is focused on EX_mnl_e.

niche_reactions.py is a statistics script rather than a COBRA one. It merges reaction presence/absence with isolation-source metadata, normalizes niche labels, runs Fisher exact tests for reaction enrichment by niche, and writes significant reaction-niche associations.