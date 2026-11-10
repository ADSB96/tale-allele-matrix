# TALE allele matrix pipeline

Scripts to predict TAL effectors (TALEs) in *Xanthomonas* genomes with
AnnoTALE, group them into classes, and build a strain × TALE-class **allele
matrix** from their RVD (repeat variable diresidue) sequences.

Within each TALE class (e.g. `TalBF`, `TalCD`), every distinct RVD sequence is
an allele. Alleles are numbered by how many strains carry them: allele 1 is the
most common.

## Repository contents

| File | Description |
|---|---|
| `annotale_pipeline.sh` | Step 1: shell commands to run AnnoTALE on a folder of genomes, assign TALEs to classes, and merge per-genome results into `AFXoc_TAL_RVD.fa` and `AFXoc_TAL_Names.tsv`. |
| `make_allelemat.R` | Step 2: R script that reads the two merged files and builds the allele matrix plus a per-TALE RVD table. |

## Workflow overview

```
Genomes/*.fasta
      │  annotale_pipeline.sh  (AnnoTALE predict → analyze → class builder → assign → rename)
      ▼
AFXoc_TAL_RVD.fa  +  AFXoc_TAL_Names.tsv
      │  make_allelemat.R
      ▼
RVD_table.xlsx  +  allelemat_new.xlsx
```

## Step 1: Predict and classify TALEs (`annotale_pipeline.sh`)

### Requirements
- Java
- [AnnoTALE](http://www.jstacs.de/index.php/AnnoTALE) command-line tool
  (`AnnoTALEcli-1.5.jar`) in the working folder

### Setup
- Put all genome assemblies (`.fasta`) in one folder called `Genomes`.
- Replace parentheses in genome file names with `_`.
- Run the commands from a folder **next to** `Genomes` (the loops read
  `../Genomes/*.fasta`). AnnoTALE creates one output folder per genome there.

### Steps
1. **predict**: find TALEs in each genome.
2. **analyze**: split each TALE into parts and RVDs.
3. **loadAndView**: interactive step. Build TALE classes and save them as
   `CLASS/Class_builder_download.xml`.
4. **assign**: assign each genome's TALEs to the classes from step 3.
5. **rename**: rename TALE predictions (GenBank, GFF3 and RVD FASTA) using the
   assigned class names.
6. **Merge files**: add the strain name to every sequence header and combine
   all genomes into:
   - `AFXoc_TAL_RVD.fa`: RVD sequences of all TALEs
   - `AFXoc_TAL_Names.tsv`: TALE ID → class name key
   - `AFXoc_TAL_AA.fa`, `AFXoc_TAL_DNA.fa`: protein and DNA sequences
   - `AFXoc_TAL_Locations.tsv`: TALE loci (GFF3)

The script is a record of the commands that were run, block by block. Run the
blocks one at a time instead of executing the whole file at once, because step 3
is interactive.

## Step 2: Build the allele matrix (`make_allelemat.R`)

### Requirements
- R ≥ 4.1
- Packages: `dplyr`, `tidyr` (≥ 1.2), `stringr`, `tibble`
- Optional: `writexl` to write `.xlsx` files (without it, `.csv` files are
  written instead)

```r
install.packages(c("dplyr", "tidyr", "stringr", "tibble", "writexl"))
```

### Input files

**`AFXoc_TAL_RVD.fa`**: one record per TALE. The header starts with
`<strain>|<TALE id>`; a `(Pseudo)` tag and bracketed coordinates may follow.
RVDs are separated by `-`.

```
>BAI119|tempTALE1 [177042-180735:1]
NI-HD-NI-N*-NI-HG-NI-NI-N*-NG-HD-NN-NI-NG-HD-NN-NS-NN-NG-NG
```

**`AFXoc_TAL_Names.tsv`**: tab-separated, no header. Column 1 is the TALE ID
exactly as in the FASTA header (without coordinates). Column 2 is the TALE name.
The class is the name without its trailing digits (`TalBU18` → `TalBU`).

```
BAI119|tempTALE1	TalBU18
BAI119|tempTALE18 (Pseudo)	TalCC8
```

### Configure and run
Edit the **USER SETTINGS** block at the top of `make_allelemat.R`:

| Setting | Purpose |
|---|---|
| `fasta_file` | Path to the RVD FASTA file |
| `names_file` | Path to the TALE name key (TSV) |
| `rvd_table_file` | Output: per-TALE table |
| `allelemat_file` | Output: allele matrix |
| `exclude_strains` | Strains to leave out (e.g. an outgroup) |
| `strain_rename` | Fix strain-name spellings, as `c("as_in_file" = "correct")` |
| `strain_order` | Row order; also breaks frequency ties. `NULL` = alphabetical |
| `class_order` | Column order. `NULL` = most widespread class first |

Then run it from the folder that contains the input files:

```bash
Rscript make_allelemat.R
```

### Allele numbering rules
- Each distinct RVD string within a class is one allele. Pseudogenes are
  included.
- Frequency is the number of **distinct strains** carrying the allele.
- Ties go to the allele found in the strain that comes earliest in
  `strain_order`, and then to the one that comes first in the FASTA file.
- A strain with two **different** alleles of a class gets `a/b` (smallest
  first, e.g. `1/3`). Two identical copies get a single number.
- `NA` means the strain has no TALE of that class.

### Outputs

**`allelemat_new.xlsx`**: one row per strain, one column per TALE class.

| strain | TalBF | TalBI | TalBO | … |
|---|---|---|---|---|
| BLS256 | 4 | 4 | 4 | … |
| CFBP7331 | 1/3 | 3 | 3 | … |
| BAI12 | 1 | 2/5 | 1 | … |

**`RVD_table.xlsx`**, with two sheets:
- `RVD_table`: one row per TALE. Columns are strain, TALE ID, name, class,
  allele name (e.g. `TalBF_Allele1`), allele number, number of strains with
  that allele, pseudogene flag, RVD count, and then `RVD1`…`RVDn`.
- `Allele_summary`: one row per allele, giving its class, number, frequency
  and full RVD string.

The script also reports IDs that are in only one input file, classes that are
listed in `class_order` but missing from the data (or the reverse), and summary
counts.
