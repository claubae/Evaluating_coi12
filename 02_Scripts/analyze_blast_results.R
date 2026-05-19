# Este script procesa los resultados de BLAST de la base de datos 12.
# Para cada archivo de salida de BLAST:
#   1. Aplica un filtro de calidad base compartido 
#   2. Guarda los top-5 hits (para ek boxplot de pident)
#   3. Aplica LCA sobre esos top 5 (si existen top 5) para comparar a asignaciones taxonómicas.


#Configurar entorno (librerías + parámetros)
library(tidyverse)
library(insect)
library(here)

source(here("betteringcoi/R/calculate_lca.R")) 

# Ejemplo: Rscript analyze_blast_results.R path/to/output_123_XX001_blocked_456.txt 197 -> así no hace falta ajustar este parámetro cada vez que cambiamos de db.
params           <- commandArgs(TRUE)
blast_output     <- params[1]
min_qlen         <- as.numeric(params[2])   # 197 para base 12, 297 para base 123

# Filtros
min_length_prop  <- 0.85
pident_abs_limit <- 80
max_evalue       <- 1e-6
min_pident_rel   <- 0.99
n_top            <- 5

file_id <- tools::file_path_sans_ext(basename(blast_output))
message("Procesando: ", blast_output)

#Extraemos el q_taxID, el b_taxID y el accession del propio archivo. 
metadata <- tibble(file = blast_output) |>
  extract(file,
          into    = c("query_taxid", "accession", "blocked_taxid"),
          regex   = "output_(\\d+)_([A-Z]+[0-9]+)_blocked_(\\d+)",
          convert = TRUE)

q_taxid   <- metadata$query_taxid
accession <- metadata$accession
b_taxid   <- metadata$blocked_taxid

# Cargar taxonomía verdadera.
taxonomy_db   <- read_csv(here("taxonomy_now.csv"), show_col_types = FALSE)
taxonomy_cols <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
assigned_cols <- paste0("assigned_", taxonomy_cols)

# Rango taxonómico del taxID bloqueado
blocked_info <- taxonomy_db |>
  filter(taxID == b_taxid) |>
  select(blocked_rank = rank)

if (nrow(blocked_info) == 0) {
  blocked_info <- tibble(blocked_rank = "unknown")
}

#Linaje de los staxIDs (obtenidos con taxonkit)
all_lineages <- read_lines(here("03_blast_results/12_good/lineages_12.txt")) |>
  tibble(raw = _) |>
  extract(raw, into = c("taxID", "taxonomy"), regex = "^(\\d+)\\s+(.*)$") |>
  separate(taxonomy, into = taxonomy_cols, sep = ";", fill = "right") |>
  mutate(across(everything(), str_trim)) |>
  rename_with(~ paste0("assigned_", .x), .cols = taxonomy_cols) |>
  mutate(taxID = as.numeric(taxID)) |>
  distinct(taxID, .keep_all = TRUE)

# Leer BLASTs (10 columnas)
blast_raw <- read_delim(
  blast_output,
  delim      = "\t",
  col_names  = c("qseqid", "sseqid", "qlen", "slen", "pident",
                 "mismatch", "evalue", "bitscore", "length", "staxids"),
  col_types  = cols(.default = col_character()),
  show_col_types = FALSE
) |>
  separate_rows(staxids, sep = ";") |>        # un staxID por fila
  mutate(
    qseqid = as.character(q_taxid),
    across(c(qlen, slen, pident, length, evalue, bitscore, staxids), as.numeric)
  ) |>
  add_column(blocked_rank = blocked_info$blocked_rank, .before = 3) |> 
  add_column(accession = metadata$accession, .before = 2)

# Unir con linajes
blast_with_taxa <- blast_raw |>
  inner_join(all_lineages, by = c("staxids" = "taxID"))

message("Resultados totales iniciales: ", nrow(blast_with_taxa))

# Calidad base (antes del LCA)
quality_filtered <- blast_with_taxa |>
  drop_na(pident, evalue, length, qlen) |>
  filter(
    length >= (min_length_prop * qlen),
    pident >= pident_abs_limit,
    evalue <= max_evalue,
    qlen   >= min_qlen
  )

message("Nº resultados tras filtro de calidad: ", nrow(quality_filtered))

# Calcular precisión taxonómica.
calculate_taxonomic_precision <- function(row, tax_levels) {
  real_precision <- "Worse"
  for (level in tax_levels) {
    t_v <- as.character(row[[paste0("true_", level)]])
    a_v <- as.character(row[[paste0("assigned_", level)]])
    if (length(t_v) == 0 || length(a_v) == 0) break
    if (is.na(t_v) || is.na(a_v)) break
    if (t_v == a_v) {
      real_precision <- level
    } else {
      break
    }
  }
  return(real_precision)
}

# Verdadero linaje del query
true_lineage_data <- get_lineage(q_taxid, taxonomy_db) |>
  enframe() |>
  filter(name %in% taxonomy_cols) |>
  pivot_wider(names_from = name, values_from = value) |>
  add_column(q_taxid = q_taxid, .before = 1) |>
  rename_with(~ paste0("true_", .x), .cols = -q_taxid)

if (nrow(true_lineage_data) == 0) {
  stop("El q_taxid ", q_taxid, " no pudo ser resuelto a un linaje completo.")
}

message("Analizando resultados del taxID: ", q_taxid,
        " Rango Bloqueado: ", blocked_info$blocked_rank)

write_csv(
  true_lineage_data,
  here("03_blast_results/12_good/true_lineages_all_v2", paste0("true_lineage_", file_id, ".csv"))
)

if (nrow(quality_filtered) == 0) {
  top5_hits <- tibble(
    qseqid       = as.character(q_taxid),
    blocked_rank = blocked_info$blocked_rank,
    staxids      = NA_real_,
    pident       = NA_real_,
    length       = NA_real_,
    !!!setNames(rep(NA_character_, length(taxonomy_cols)), assigned_cols)
  )
  message("No hits (top5) for ", q_taxid)
} else {
  top5_hits <- quality_filtered |>
    filter(pident >= (min_pident_rel * max(pident, na.rm = TRUE))) |>
    slice_max(order_by = pident, n = n_top, with_ties = TRUE)
}

top5_combined <- top5_hits |>
  mutate(qseqid = as.character(qseqid)) |>
  left_join(
    true_lineage_data |> mutate(q_taxid = as.character(q_taxid)),
    by = c("qseqid" = "q_taxid")
  ) |>
  mutate(across(everything(), as.character))

top5_combined$precision <- pmap_chr(
  top5_combined,
  function(...) calculate_taxonomic_precision(list(...), taxonomy_cols)
)

write_csv(
  top5_combined,
  here("03_blast_results/12_good/individual_results_all_v2", paste0("res_", file_id, ".csv"))
)
message("Nº de top 5 hits guardados: ", nrow(top5_combined))

#Aquí obtenemos resultados para la gráfica comparativa del pident

#Ahora calculamos el LCA de los mejores hits.
lca_result <- calculate_lca(top5_hits)

if (nrow(lca_result) == 0) {
  final_result <- tibble(
    qseqid       = as.character(q_taxid),
    blocked_rank = blocked_info$blocked_rank,
    staxids      = NA_real_,
    pident       = NA_real_,
    length       = NA_real_,
    !!!setNames(rep("No_Hit", length(taxonomy_cols)), assigned_cols)
  )
  message("No hits (LCA) for ", q_taxid)
} else {
  final_result <- lca_result
}

#Se vuelve a calcular la precisión taxonómica (ya que puede cambiar tras LCA)
combined_data <- final_result |>
  mutate(qseqid = as.character(qseqid)) |>
  left_join(
    true_lineage_data |> mutate(q_taxid = as.character(q_taxid)),
    by = c("qseqid" = "q_taxid")
  ) |>
  mutate(across(everything(), as.character))

combined_data$precision <- pmap_chr(
  combined_data,
  function(...) {
    calculate_taxonomic_precision(list(...), taxonomy_cols)
  }
)

write_csv(
  combined_data,
  here("03_blast_results/12_good/individual_results_v2", paste0("res_", file_id, ".csv"))
)

message("Cálculo taxonomía acabado para: ", q_taxid)

#Con esto se hace comparativa entre precisión taxonómica 12 y precisión taxonómica 123.

#out_folder=/home/cbaeyens/03_blast_results/blast_parsed_12_v2
## blast_folder=/home/cbaeyens/03_blast_results/12_good
#sbatch --array=1-909 --export=BLAST_FOLDER="$blast_folder" /home/cbaeyens/01_Scripts/Step2_taxonomy_array.sh 
