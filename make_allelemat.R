#!/usr/bin/env Rscript
# =============================================================================
# make_allelemat.R
#
# Build a TALE allele matrix (strains x TALE classes) from RVD sequences.
#
# Within each TALE class (e.g. TalBF, TalCD), every distinct RVD sequence is
# treated as one allele. Alleles are numbered by how many strains carry them:
# allele 1 is found in the most strains, allele 2 in the next most, and so on.
# The final matrix has one row per strain and one column per TALE class; each
# cell holds the allele number that strain carries, "a/b" if the strain carries
# two different alleles of that class, or "NA" if the class is absent.
#
# USAGE
#   1. Edit the "USER SETTINGS" section below (input/output file names and,
#      optionally, strain/class ordering).
#   2. Run from the folder containing the input files:
#        Rscript make_allelemat.R
#
# REQUIREMENTS
#   R >= 4.1 (uses the native |> pipe)
#   Packages: dplyr, tidyr (>= 1.2, for separate_wider_delim), stringr, tibble
#   Optional: writexl  -> writes .xlsx output. If it is not installed, the
#             script writes .csv files instead.
#             Install with: install.packages(c("dplyr","tidyr","stringr","tibble","writexl"))
#
# NUMBERING RULES
#   - Each distinct RVD string within a class = one allele. Pseudogenes are
#     included and treated like any other TALE.
#   - Allele frequency = number of DISTINCT strains carrying the allele
#     (a strain with two identical copies counts once).
#   - Ties in frequency are broken by the first strain (in `strain_order`)
#     that carries the allele, then by position in the FASTA file.
#   - A strain with two identical copies of an allele gets a single number
#     (e.g. "1", not "1/1"); two different alleles are written smallest first
#     (e.g. "1/3").
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

# =============================================================================
# USER SETTINGS -- edit this section
# =============================================================================

# -----------------------------------------------------------------------------
# INPUT FILE 1: RVD sequences in FASTA format
# -----------------------------------------------------------------------------
# One record per TALE. The header must start with "<strain>|<TALE id>".
# Anything in square brackets at the end of the header (e.g. genome
# coordinates) is ignored. A "(Pseudo)" tag may be present and is kept as part
# of the ID, so it must match the ID in the names file exactly.
# The sequence is the RVDs joined by "-" (may span multiple lines).
#
# Example:
#   >BAI119|tempTALE1 [177042-180735:1]
#   NI-HD-NI-N*-NI-HG-NI-NI-N*-NG-HD-NN-NI-NG-HD-NN-NS-NN-NG-NG
#   >BAI119|tempTALE18 (Pseudo) [1582342-1583725:-1]
#   NN-HD-...
fasta_file <- "AFXoc_TAL_RVD.fa"

# -----------------------------------------------------------------------------
# INPUT FILE 2: TALE name key, tab-separated, NO header row
# -----------------------------------------------------------------------------
# Column 1: TALE ID, exactly as it appears in the FASTA header (without the
#           bracketed coordinates), e.g. "BAI119|tempTALE1" or
#           "BAI119|tempTALE18 (Pseudo)".
# Column 2: TALE name, e.g. "TalBU18". The class is the name with trailing
#           digits removed ("TalBU18" -> class "TalBU").
#
# Example:
#   BAI119|tempTALE1<TAB>TalBU18
#   BAI119|tempTALE2<TAB>TalCD9
#
# TALEs present in the FASTA but missing from this file (or vice versa) are
# reported as a warning and left out of the output.
names_file <- "AFXoc_TAL_Names.tsv"

# -----------------------------------------------------------------------------
# OUTPUT FILES (written to the current working directory)
# -----------------------------------------------------------------------------
# Intermediate per-TALE table. Two sheets:
#   "RVD_table"      - one row per TALE: strain, ID, name, class, allele name
#                      and number, number of strains with that allele,
#                      pseudogene flag, RVD count, then one RVD per column.
#   "Allele_summary" - one row per allele: class, allele number, number of
#                      strains carrying it, and its full RVD string.
rvd_table_file <- "RVD_table.xlsx"

# Final allele matrix: strains as rows, TALE classes as columns.
allelemat_file <- "allelemat_new.xlsx"

# -----------------------------------------------------------------------------
# STRAIN CLEAN-UP
# -----------------------------------------------------------------------------
# Strains to leave out of the analysis entirely (they also do not count toward
# allele frequencies). Use character(0) to keep all strains.
exclude_strains <- c("BAI3")

# Fix strain-name spellings found in the input files: c("as_in_file" = "correct").
# Use character(0) if no renaming is needed.
strain_rename <- c("CBFP7331" = "CFBP7331",
                   "CBFP7341" = "CFBP7341")

# -----------------------------------------------------------------------------
# ROW ORDER (strains)
# -----------------------------------------------------------------------------
# Order of strains (rows) in the output, using the corrected names from above.
# This order is also used to break ties in allele frequency.
# Strains in the data but not listed here are appended at the end (with a
# warning). Set to NULL to order strains alphabetically.
strain_order <- c(
  "BLS256", "RS105", "BXOR1", "MD-P1-1B", "CAAWA_1a", "Xoc010", "Ug_14",
  "MB-P7-3-2", "MB-P38-1-1", "Bur7", "Bur2", "Bur6", "MB-P13-1-1", "Xoc007",
  "TANZ-P11-2L", "MD-P13-1-7", "PSMORO-1", "MD-Ivory-1-B", "Selingue50",
  "Selingue51", "BAI6", "S23-1-10", "S52-1-6", "CFBP7342", "BAI12", "MAI35",
  "CFBP7331", "MAI23", "CFBP7341", "BAI375", "BAI7", "MAI25", "MAI7",
  "BAI119", "MAI42"
)

# -----------------------------------------------------------------------------
# COLUMN ORDER (TALE classes)
# -----------------------------------------------------------------------------
# Order of TALE classes (columns) in the output. Classes listed here but not
# found in the data are dropped; classes in the data but not listed here are
# appended at the end in alphabetical order (both are reported as messages).
# Set to NULL to order classes by the number of strains carrying them
# (most common first).
class_order <- c(
  "TalBF", "TalBI", "TalBO", "TalBD", "TalBU", "TalBC", "TalAT", "TalAY",
  "TalBW", "TalBE", "TalAV", "TalBL", "TalAX", "TalCP", "TalCD", "TalAZ",
  "TalBN", "TalCE", "TalCG", "TalAD", "TalCL", "TalBR", "TalCO", "TalBV",
  "TalHK", "TalID", "TalCB", "TalCJ", "TalCF", "TalCH", "TalCI", "TalCC",
  "TalCR", "TalBB", "TalAK", "TalCN", "TalCQ", "TalJA", "TalJE", "TalHL",
  "TalCM", "TalBG", "TalBY", "TalAJ", "TalAU", "TalJU", "TalKX", "TalJR",
  "TalKT", "TalKS", "TalKY"
)

# =============================================================================
# END OF USER SETTINGS -- no edits needed below this line
# =============================================================================


# ---- 1. Read the FASTA file --------------------------------------------------
# Returns one row per record: the full header, the TALE ID (header minus the
# bracketed coordinates) and the RVD string (multi-line sequences are joined).
read_rvd_fasta <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  lines  <- lines[nzchar(trimws(lines))]          # drop blank lines
  is_hdr <- startsWith(lines, ">")
  grp    <- cumsum(is_hdr)                        # record number for each line
  hdr    <- sub("^>", "", lines[is_hdr])
  seqs   <- vapply(split(lines[!is_hdr], grp[!is_hdr]), paste, "", collapse = "")
  tibble(
    header = hdr,
    id     = str_trim(str_remove(hdr, "\\s*\\[.*\\]\\s*$")),
    rvd    = unname(seqs[as.character(seq_along(hdr))])
  )
}

fa <- read_rvd_fasta(fasta_file)

# ---- 2. Read the TALE name key -----------------------------------------------
tsv <- read.delim(names_file, header = FALSE, sep = "\t", quote = "",
                  col.names = c("id", "name"), stringsAsFactors = FALSE) |>
  as_tibble() |>
  mutate(id = str_trim(id), name = str_trim(name))

# Report IDs that are present in only one of the two input files.
missing_in_tsv <- setdiff(fa$id, tsv$id)
missing_in_fa  <- setdiff(tsv$id, fa$id)
if (length(missing_in_tsv)) warning("In FASTA but not in names file (skipped): ",
                                    paste(missing_in_tsv, collapse = ", "))
if (length(missing_in_fa))  warning("In names file but not in FASTA (skipped): ",
                                    paste(missing_in_fa,  collapse = ", "))

# ---- 3. Join and annotate each TALE ------------------------------------------
tales <- inner_join(fa, tsv, by = "id") |>
  mutate(
    strain   = str_extract(id, "^[^|]+"),          # text before the "|"
    strain   = recode(strain, !!!strain_rename),   # apply spelling fixes
    pseudo   = str_detect(id, "\\(Pseudo\\)"),
    class    = str_remove(name, "\\d+$"),          # TalBU18 -> TalBU
    n_rvd    = str_count(rvd, "-") + 1,
    file_pos = row_number()                        # used as a final tie-breaker
  ) |>
  filter(!strain %in% exclude_strains)

# Resolve strain order: default alphabetical; append any unlisted strains.
if (is.null(strain_order)) strain_order <- sort(unique(tales$strain))
unknown_strains <- setdiff(unique(tales$strain), strain_order)
if (length(unknown_strains)) {
  warning("Strains not in strain_order (appended at end): ",
          paste(unknown_strains, collapse = ", "))
  strain_order <- c(strain_order, unknown_strains)
}
tales <- tales |> mutate(strain_rank = match(strain, strain_order))

# ---- 4. Define and number alleles --------------------------------------------
# One row per (class, RVD string). Sort by frequency (desc), then by the first
# strain carrying it, then by FASTA position, and number within each class.
alleles <- tales |>
  group_by(class, rvd) |>
  summarise(
    n_strains  = n_distinct(strain),
    first_rank = min(strain_rank),
    first_pos  = min(file_pos[strain_rank == min(strain_rank)]),
    .groups = "drop"
  ) |>
  arrange(class, desc(n_strains), first_rank, first_pos) |>
  group_by(class) |>
  mutate(allele = row_number()) |>
  ungroup()

tales <- tales |>
  left_join(select(alleles, class, rvd, allele, n_strains), by = c("class", "rvd")) |>
  mutate(allele_name = paste0(class, "_Allele", allele))

# ---- 5. Resolve column (class) order -----------------------------------------
cols_present <- unique(tales$class)
if (is.null(class_order)) {
  class_order <- tales |>
    distinct(class, strain) |>
    count(class, name = "n") |>
    arrange(desc(n), class) |>
    pull(class)
}
col_order <- c(intersect(class_order, cols_present),
               sort(setdiff(cols_present, class_order)))

dropped_cols <- setdiff(class_order, cols_present)
if (length(dropped_cols)) message("Classes in class_order not found in data (dropped): ",
                                  paste(dropped_cols, collapse = ", "))
added_cols <- setdiff(col_order, class_order)
if (length(added_cols)) message("Classes in data not in class_order (appended): ",
                                paste(added_cols, collapse = ", "))

# ---- 6. Intermediate per-TALE table ------------------------------------------
# Splits the RVD string into one column per repeat (RVD1, RVD2, ...).
rvd_cols <- paste0("RVD", seq_len(max(tales$n_rvd)))

rvd_table <- tales |>
  arrange(match(class, col_order), allele, strain_rank) |>
  select(Strain = strain, TALE_ID = id, Name = name, Class = class,
         AlleleName = allele_name, Allele = allele, StrainsWithAllele = n_strains,
         Pseudo = pseudo, nRVD = n_rvd, rvd) |>
  separate_wider_delim(rvd, delim = "-", names = rvd_cols, too_few = "align_start")

allele_summary <- alleles |>
  arrange(match(class, col_order), allele) |>
  select(Class = class, Allele = allele, StrainsWithAllele = n_strains, RVD = rvd)

# ---- 7. Allele matrix --------------------------------------------------------
# Collapse to one value per strain x class ("1", or "1/3" for two alleles),
# then pivot to wide format with strains as rows.
cells <- tales |>
  distinct(strain, class, allele) |>
  arrange(strain, class, allele) |>
  group_by(strain, class) |>
  summarise(value = paste(allele, collapse = "/"), .groups = "drop")

allelemat <- cells |>
  pivot_wider(names_from = class, values_from = value) |>
  right_join(tibble(strain = strain_order), by = "strain") |>
  arrange(match(strain, strain_order)) |>
  select(strain, all_of(col_order)) |>
  mutate(across(-strain, ~ replace_na(.x, "NA")))

# ---- 8. Write outputs --------------------------------------------------------
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(RVD_table = rvd_table, Allele_summary = allele_summary),
                      rvd_table_file)
  writexl::write_xlsx(list(allelemat = allelemat), allelemat_file)
  message("Wrote ", rvd_table_file, " and ", allelemat_file)
} else {
  # Fallback: same tables as CSV, named after the .xlsx settings above.
  base_rvd <- tools::file_path_sans_ext(rvd_table_file)
  base_mat <- tools::file_path_sans_ext(allelemat_file)
  write.csv(rvd_table,      paste0(base_rvd, ".csv"),         row.names = FALSE, na = "")
  write.csv(allele_summary, paste0(base_rvd, "_summary.csv"), row.names = FALSE)
  write.csv(allelemat,      paste0(base_mat, ".csv"),         row.names = FALSE)
  message("writexl not installed; wrote CSV files instead")
}

message(sprintf(
  "%d TALEs, %d strains, %d classes, %d alleles; %d strain/class cells with 2+ alleles",
  nrow(tales), n_distinct(tales$strain), length(col_order), nrow(alleles),
  sum(str_detect(cells$value, "/"))
))
