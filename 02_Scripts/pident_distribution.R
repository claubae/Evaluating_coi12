#Este script filtra los resultados de BLAST y los categoriza entre "precisión taxonómica es inmediatamente superior al rango bloqueado" o "Worse"
#Se crea un df nuevo con la mediana, por cada query, del pident cuando precisión es rango inm. superior y mediana de cuando es Worse.

#Configurar entorno (librerías + parámetros)
library(tidyverse)
library(insect)
library(here)
params           <- commandArgs(TRUE)
blast_output     <- params[1]
min_qlen         <- params[2] # 197 para base 12, 297 para base 123 -> se especifica en llamamiento script

# Filtros 
min_length_prop  <- 0.85
pident_abs_limit <- 80
max_evalue       <- 1e-6


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

# Determinar rango taxonómico del taxID bloqueado
blocked_info <- taxonomy_db |>
  filter(taxID == b_taxid) |>
  select(blocked_rank = rank)

if (nrow(blocked_info) == 0) {
  blocked_info <- tibble(blocked_rank = "no_block")
} #los que no tienen ningún rango bloqueado se llaman -> ..._blocked_0000.tsv.

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
  separate_rows(staxids, sep = ";") |>      
  mutate(
    qseqid = as.character(q_taxid), 
    across(c(qlen, slen, pident, length, evalue, bitscore, staxids), as.numeric)
  ) |>
  add_column(blocked_rank = blocked_info$blocked_rank, .before = 3) |> 
  add_column(accession = metadata$accession, .before = 2)

# Unir con linajes asignados (determinados con taxonkit)
blast_with_taxa <- blast_raw |>
  inner_join(all_lineages, by = c("staxids" = "taxID"))

message("Resultados totales iniciales: ", nrow(blast_with_taxa))

# Calidad base
quality_filtered <- blast_with_taxa |>
  drop_na(pident, evalue, length, qlen) |>
  filter(
    length >= (min_length_prop * qlen),
    pident >= pident_abs_limit,
    evalue <= max_evalue,
    qlen   >= min_qlen
  )

message("Nº resultados tras filtro de calidad: ", nrow(quality_filtered))

# Definir función para calcular precisión taxonómica.
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

if (nrow(quality_filtered) == 0) {
  pident_filter <- tibble(
    qseqid      = as.character(q_taxid),
    blocked_rank = blocked_info$blocked_rank,
    staxids      = NA_real_,
    pident       = NA_real_,
    length       = NA_real_,
    !!!setNames(rep(NA_character_, length(taxonomy_cols)), assigned_cols)
  )
  message("No hits for ", q_taxid)
} else {
  pident_filter <- quality_filtered 
}

pident_combined <- pident_filter |>
  mutate(qseqid = as.character(qseqid)) |>
  left_join(
    true_lineage_data |> mutate(q_taxid = as.character(q_taxid)),
    by = c("qseqid" = "q_taxid")
  ) |>
  mutate(across(everything(), as.character)) |> 
  mutate(unique_id = paste0(q_taxid, "_", accession, "_block_", blocked_rank), .before = 2)
#Ahora calculamos precisión taxonómica
pident_combined$precision <- pmap_chr(
  pident_combined,
  function(...) calculate_taxonomic_precision(list(...), taxonomy_cols)
)

#Definimos los rangos taxonómicos para poder categorizarlo en "Worse" o precisión en rango inmediatamente superior al bloqueado:
blocked_to_next_level <- c(
  "no_block" = "species",
  "species" = "genus",
  "genus"   = "family",
  "family"  = "order",
  "order"   = "class",
  "class"   = "phylum"
)

precision_levels <- c("species", "genus", "family", "order", "class", "phylum", "kingdom", "Worse")

pident_combined <- pident_combined |> 
  mutate(
    pident       = as.numeric(pident),
    blocked_rank = factor(blocked_rank, levels = names(blocked_to_next_level)),
    next1        = unname(blocked_to_next_level[as.character(blocked_rank)]),
    precision_group = if_else(precision == next1, next1, "Worse"),
    precision_group = factor(precision_group, levels = c("species", "genus", "family", "order", "class", "phylum", "Worse"))
  ) |> 
  add_column(database="12", .before = 2)

write_csv(
  pident_combined,
  here("03_blast_results/12_good/pident_calculations", paste0("pident_", file_id, ".csv"))
)

#Ahora podemos pedir que nos una todos los archivos en uno solo:

pident_path <- here("03_blast_results/12_good/pident_calculations")

pident_files <- list.files(path = pident_path, pattern = "\\.csv$", full.names = TRUE)

all_pident <- pident_files |> 
  set_names(tools::file_path_sans_ext(basename(pident_files))) |> 
  map_dfr(function(f) read_csv(f, col_types = cols(.default = "c"), show_col_types = FALSE))

sample_size <- nrow(all_pident)

#Una vez que los tenemos todas las obs. juntas pasamos pident -> prop para poder hacer la regresión beta mixta -> ordbeta para que entre interval cerrado [0;1]
#Ya que tenemos valores de pident=100%-> ordbeta (Kubinec, 2022)
pident_mod <- all_pident |>
  mutate(
    pident = as.numeric(pident),
    #Leer documentación glmmTMB para ordbeta (dentro de betafamily) -> transformar intervalo [a, b] a [0,1] -> (y-a/b-a), siendo nuestro intervalo [80;100]
    y_norm = ((pident-80)/(100-80)),
    precision_group=fct_relevel(precision_group, "species", "genus", "family", "order", "class", "phylum", "Worse"),
    status= if_else(precision == next1, "Correcto", "Incorrecto"))
saveRDS(pident_mod, file = here("03_blast_results/12_good/pident_calculations/all_pident_together.rds"))


## out_folder=/home/cbaeyens/03_blast_results/pident_123_parsed
## blast_folder=/home/cbaeyens/03_blast_results/123_good
## sbatch --array=1-909 --export=BLAST_FOLDER="$blast_folder" /home/cbaeyens/01_Scripts/Step2_taxonomy_array.sh 