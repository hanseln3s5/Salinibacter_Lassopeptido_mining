# Genome mining of *Salinibacter ruber*: prioritisation of biosynthetic gene clusters and *in silico* structural modelling of lasso peptides

Analysis code for the master's thesis *Minería genómica de Salinibacter ruber:
priorización de clústeres biosintéticos y modelado estructural in silico de
lassopéptidos con potencial farmacológico* (Máster Universitario en
Bioinformática y Bioestadística, UOC–UB).

Author: Hansel B. López-Chevalier
Licence: MIT

---

## What this repository is for

Every script here exists to close a specific gap that the standard tools leave
open. The reasoning behind each design decision is documented in the module
docstring of the script itself, not only here, so that a reviewer reading the
code never has to guess why a threshold was chosen.

The repository does **not** reimplement antiSMASH, RODEO, LassoPred, LassoProt
or HADDOCK. Those are used as published. What is provided here is the
connective tissue between them, plus the checks that none of them perform.

---

## Pipeline

Numbered scripts run in order. External tools are marked.

| Step | Tool | Type | Input → Output |
|------|------|------|----------------|
| 0 | antiSMASH | external | genomes → cluster regions (`.gbk`) |
| 1 | `01_mine_bgcs.py` | this repo | Prokka output → marker-domain survey, assembly quality, false-negative flags |
| 2 | `02_find_precursors.py` | this repo | `.gbk` → ranked precursor candidates, core peptides |
| 2b | RODEO, BAGEL4 | external (web) | `.gbk` → confirmed precursor |
| 3 | LassoPred | external (web) | core sequence → `.pdb` |
| 4 | `03_validate_lasso_fold.py` | this repo | `.pdb` → topology verdict + **`.pdb` carrying the LINK record** |
| 4b | LassoProt | external (web) | LINK-bearing `.pdb` → independent classification |
| 5 | `04_characterize_peptide.py` | this repo | `.pdb` → physicochemistry, surface, geometry |
| 5b | PDB2PQR → APBS → PyMOL | external | `.pdb` → electrostatic surface |
| 6 | CHARMM-GUI | external (web) | LINK-bearing `.pdb` → MD topology with the crosslink declared |
| 6b | GROMACS or OpenMM | external | → trajectory |
| 7 | `05_analyze_md_trajectory.py` | this repo | trajectory → RMSD, RMSF, **per-frame threading** |
| 8 | `06_prepare_receptor.py` | this repo | reference complex → receptor, reference peptide, contact residues |
| 9 | `07_interface_residues.py` | this repo | complex → HADDOCK active residues, with accessibility |
| 10 | `08_prepare_haddock_ligand.py` | this repo | model(s) → ensemble `.pdb` + isopeptide restraint `.tbl` |
| 11 | HADDOCK | external (web) | → docked poses |

`aux_osmoadaptation_profile.py` sits outside the pipeline. It classifies the
osmoadaptation strategy of a genome from its proteome, and provides the
physiological context in which a strongly acidic core peptide is interpreted.

### Two ordering constraints that matter

**`03` must run before `04` and before any MD setup.** It writes the PDB `LINK`
record for the isopeptide bond. LassoProt, CHARMM-GUI and MD topology builders
do not infer non-canonical covalent bonds from geometry; without that record
the macrolactam is treated as absent, LassoProt reads the chain as linear, and
an MD topology built from the file lets the ring fall open during the
simulation. **Use the `_LINK.pdb` output of `03` for everything downstream, not
the raw predictor output.**

**`02` must run before `03`.** No core peptide, no structure to validate.

---

## Scripts

### `01_mine_bgcs.py`
antiSMASH-independent survey of lasso peptide marker domains across a genome
set, plus assembly-quality statistics. Answers whether negative genomes are
truly negative or whether the cluster was missed because of fragmentation or a
contig edge. `PF13471` is treated as the discriminating marker; `PF00733`
counts only when co-localised with it, because housekeeping asparagine
synthetase is present in nearly every genome.

Requires HMMER on `PATH`. Fetch the HMMs with `--download-hmm hmms/`.

### `02_find_precursors.py`
Recovers lasso peptide precursor genes that general-purpose gene callers miss.
Prodigal imposes a 90 nt minimum ORF length; lasso precursors are 25–60
residues and often start at GTG or TTG, so they are never annotated. Performs a
six-frame scan honouring alternative start codons, then scores candidates on
genomic context, ring-forming residue, ring size, core length, tail and plug,
and leader composition.

**This is a heuristic screen, not a replacement for RODEO.** Report it as an
exploratory complementary analysis in an annex, never alongside the trained
tools. See the module docstring for the full caveat.

### `03_validate_lasso_fold.py`
Verifies that the predicted model actually is a lasso peptide: measures the
isopeptide bond and counts how many times the tail crosses the macrocycle, by
Möller–Trumbore ray–triangle intersection against a triangle fan over the ring
Cα atoms. An odd count means threaded. Writes the `LINK` record.

Exit status 0 when the model passes both criteria, 1 otherwise, so it can gate
a shell pipeline.

### `04_characterize_peptide.py`
Physicochemistry, surface and geometry. Corrects the isoelectric point, net
charge and molecular weight for the fact that the α-amino group of residue 1
and the acceptor side-chain carboxyl are consumed by the isopeptide bond —
which a general-purpose tool such as ProtParam gets wrong. Reports both the
naive and corrected values.

Detects and reports explicitly the case where the peptide has **no accessible
isoelectric point**, rather than returning a meaningless 0.00.

### `05_analyze_md_trajectory.py`
Trajectory analysis including the check no general tool performs: whether the
threading survives the simulation, recomputed frame by frame. RMSD tells you
the structure moved; it does not tell you it is still a lasso peptide.

Also monitors the isopeptide bond length and warns when it exceeds 2.0 Å, which
almost always means the bond was not declared in the simulation topology.

Requires MDAnalysis.

### `06_prepare_receptor.py`
Splits a reference complex (6N60: *E. coli* RNA polymerase with microcin J25),
discards nucleic acid, waters, ions and heteroatoms, renumbers residues so they
are unique across chains as HADDOCK requires, optionally truncates around the
peptide to fit the server's atom limit, and extracts the contact residues.
Always writes `receptor_numbering.tsv` so results can be reported in the
original numbering.

Requires numpy.

### `07_interface_residues.py`
Extracts interface residues in the numbering of **your own** file, with
optional Shrake–Rupley accessibility. Flags the residues HADDOCK would discard
as buried — which matters because the RNA polymerase secondary channel is a
tunnel and its lining residues are buried by definition, so leaving that option
enabled deletes the binding site.

### `08_prepare_haddock_ligand.py`
Translates AMBER residue naming to the standard names CNS understands, fixes
the segid, builds an ensemble PDB from several conformers, and writes the
unambiguous distance restraint that keeps the macrolactam closed during
docking. Without that restraint the ring opens during simulated annealing and
the molecule being docked is no longer a lasso peptide.

### `aux_osmoadaptation_profile.py`
Classifies salt-in versus salt-out strategy from proteome isoelectric point
distribution plus marker genes. Requires numpy; matplotlib only for `--figure`.

---

## Requirements

Most scripts use the Python 3.8+ standard library only, deliberately, so that
they run in a bare container with no scientific stack installed.

| Script | Dependencies |
|--------|--------------|
| `01`, `02`, `04`, `07`, `08` | none (standard library) |
| `03`, `06` | numpy |
| `05` | numpy, MDAnalysis |
| `aux_osmoadaptation_profile` | numpy; matplotlib for `--figure` |
| `01` additionally | HMMER (`hmmsearch`) on `PATH` |

```bash
conda create -n lasso python=3.11 numpy mdanalysis matplotlib
conda activate lasso
conda install -c bioconda hmmer
```

---

## Reproducibility notes

Record these when running, because they cannot be reconstructed afterwards:

- antiSMASH version and run date (both appear in the `.gbk` header).
- RODEO: commit hash, since the repository publishes no tagged releases.
- Pfam and HMMER versions used by `01`.
- The exact command line of every script, including non-default thresholds.
- For `02`, the **full** candidate table, not only the candidates above the
  score cutoff. Without it a reader cannot judge whether the cutoff was
  reasonable.
- For `06`, whether the receptor was truncated and at what radius.
- LassoPred and LassoProt access dates; both are web services that change
  without versioning.

Set a fixed `--locustag` when running Prokka. Prokka generates a random
eight-letter prefix per run, so the same genome annotated twice yields locus
tags that cannot be cross-referenced.

---

## Limitations to state in the thesis

1. `02_find_precursors.py` is heuristic and its weights were chosen from the
   literature, not trained. It produces candidates, not results.
2. RODEO, LassoPred and LassoESM were trained predominantly on Actinobacteria
   and Pseudomonadota. *Salinibacter* is a Bacteroidota, and therefore outside
   the training domain of the whole tool chain. A positive result is a strong
   lead; a negative does not prove absence.
3. LassoPred reproduces the ring considerably better than the tail. Any
   conclusion about surface or docking that depends on the exact tail
   orientation inherits that uncertainty.
4. Molecular dynamics tests whether a model is **stable** under a force field,
   not whether it is **correct**. An incorrect structure can be perfectly
   stable.
5. Docking always returns a best pose. Only comparison against a reference
   ligand of known activity distinguishes a plausible interaction from an
   artefact of the procedure.

---

## File name mapping

These scripts were originally written with Spanish names during development.
The mapping is recorded here so earlier notes and intermediate files can be
traced.

| Current | Original |
|---------|----------|
| `01_mine_bgcs.py` | `busca_lasso_bgc.py` |
| `02_find_precursors.py` | `busca_precursor_lasso.py` |
| `03_validate_lasso_fold.py` | `valida_lazo.py` |
| `04_characterize_peptide.py` | `caracteriza_lassopeptido.py` |
| `05_analyze_md_trajectory.py` | `analiza_md_lasso.py` |
| `06_prepare_receptor.py` | `prepara_receptor.py` |
| `07_interface_residues.py` | `interfaz_haddock.py` |
| `08_prepare_haddock_ligand.py` | `prepara_haddock.py` |
| `aux_osmoadaptation_profile.py` | `perfil_osmoadaptacion.py` |

Behaviour is unchanged: each translated script was run against its Spanish
original on the same input and the numerical outputs verified identical.

---

## Key references

- Tietz, J. I. *et al.* (2017). A new genome-mining tool redefines the lasso peptide biosynthetic landscape. *Nature Chemical Biology* **13**, 470–478.
- Barrett, S. E. & Mitchell, D. A. (2024). Advances in lasso peptide discovery, biosynthesis, and function. *Trends in Genetics* **40**, 950–968.
- Braffman, N. R. *et al.* (2019). Structural mechanism of transcription inhibition by lasso peptides microcin J25 and capistruin. *PNAS* **116**, 1273–1278.
- Dabrowski-Tumanski, P. *et al.* (2016). LassoProt: server to analyze biopolymers with lassos. *Nucleic Acids Research* **44**, W383–W389.
- Ouyang, X. *et al.* (2025). LassoPred: a tool to predict the 3D structure of lasso peptides. *JCIM*.
- Charitou, V., van Keulen, S. C. & Bonvin, A. M. J. J. (2022). Cyclization and docking protocol for cyclic peptide–protein modeling using HADDOCK2.4. *JCTC* **18**, 4027–4040.
- Pašić, L. *et al.* (2009). Metagenomic islands of hyperhalophiles: the case of *Salinibacter ruber*. *BMC Genomics* **10**, 570.
- Mongodin, E. F. *et al.* (2005). The genome of *Salinibacter ruber*. *PNAS* **102**, 18147–18152.
