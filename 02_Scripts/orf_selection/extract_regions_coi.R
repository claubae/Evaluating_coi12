#  extract_regions_coi.R

#  Extrae 3 regiones de COI de cada secuencia query y genera sus versiones
#  123 (completa) y 12 (sin 3ª base, con el marco de lectura de cada secuencia)

#  REGIONES (como en Leray et al. 2022, MIDORI2)
#   1. complete : gen COI completo (~1500 pb)        -> la secuencia entera
#   2. barcode  : región Folmer (~658 pb)            -> entre LCO1490 y HCO2198
#   3. leray    : fragmento Leray (~313 pb)          -> entre mlCOIintF y HCO2198
#   Los primers NO se incluyen en la región (igual que con cutadapt en el TFG).
#
#        LCO1490          mlCOIintF                       HCO2198
#   5' ===>|----------------===>|---------------------------|<=== 3'
#          |<---------------- barcode (~658) -------------->|
#                               |<------ leray (~313) ----->|
#   |<------------------------- complete (~1500) ------------------------->|

#  IMPORTANTE (diseño experimental)
#   - Para que la comparación sea justa, cada query debe tener LAS TRES regiones
#     (misma secuencia, misma especie). Por eso necesitas secuencias de COI
#     COMPLETO (p. ej. de mitogenomas) de Calanoida, sin recortar con cutadapt.
#   - La base de referencia es la MISMA para las 3 regiones: MIDORI sin recortar
#     (versión 123) y su versión 12 hecha con orf_selection_coi12.R.
#
#  USO: abre Evaluating_coi12.Rproj, cambia el PASO 1 y ejecuta bloque a bloque.
#  Salida (out_dir): <region>_123.fasta, <region>_12.fasta, regions_report.tsv



# -----------------------------------------------------------------------------
# (0) Librerías
# -----------------------------------------------------------------------------

  library(Biostrings)
  library(here)
  library(readr)
  library(dplyr)



# -----------------------------------------------------------------------------
# (1) PARÁMETROS
# -----------------------------------------------------------------------------
input_file    <- here("data", "calanoida_coi_completo.fasta")  # FASTA o CSV (taxid,seq,accession,nombre)
out_dir       <- here("08_orf_results", "regiones_coi")
genetic_table <- 5        # mitocondrial invertebrados
max_mismatch  <- 4        # errores permitidos en cada primer (como -e 4 en cutadapt)
min_complete  <- 1400     # longitud mínima para considerar el COI "completo"
same_queries  <- TRUE     # TRUE = solo queries con las 3 regiones (comparación justa)

# Primers (IUPAC; las inosinas "I" se escriben como N)
primers <- list(
  LCO1490   = "TNTCNACNAAYCAYAARGAYATTGG",   # jgLCO1490 (Geller et al. 2013)
  mlCOIintF = "GGWACWGGWTGAACWGTWTAYCCYCC",  # Leray et al. 2013 (el del TFG)
  HCO2198rc = "TGRTTYTTYGGNCAYCCNGARGTNTA"   # reverso complementario de jgHCO2198
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# (2) Leer secuencias (FASTA o CSV del proyecto)
# -----------------------------------------------------------------------------
read_input <- function(file) {
  if (substr(readLines(file, n = 1), 1, 1) == ">") {
    s <- readDNAStringSet(file)
  } else {
    tab <- read_csv(file, col_names = c("taxid", "seq", "accession", "scientific_name"),
                    show_col_types = FALSE)
    s <- DNAStringSet(tab$seq)
    names(s) <- paste0(tab$accession, ";", tab$scientific_name)
  }
  DNAStringSet(gsub("-", "", toupper(as.character(s))))
}

seqs <- read_input(input_file)
cat("Secuencias leídas:", length(seqs), "\n")
summary(width(seqs))


# -----------------------------------------------------------------------------
# (3) Buscar los primers en cada secuencia
# -----------------------------------------------------------------------------
# fixed = FALSE -> los códigos IUPAC (W, Y, R, N...) se interpretan.
# Devuelve la posición del PRIMER acierto (inicio y fin) o NA si no aparece.
find_primer <- function(seqs, primer, mm, from = rep(1L, length(seqs))) {
  start <- end <- rep(NA_integer_, length(seqs))
  for (i in seq_along(seqs)) {
    if (is.na(from[i]) || from[i] > width(seqs)[i]) next
    sub <- subseq(seqs[[i]], start = from[i])
    hit <- matchPattern(DNAString(primer), sub, max.mismatch = mm, fixed = FALSE)
    if (length(hit) > 0) {
      start[i] <- start(hit)[1] + from[i] - 1L
      end[i]   <- end(hit)[1]   + from[i] - 1L
    }
  }
  data.frame(start = start, end = end)
}

lco   <- find_primer(seqs, primers$LCO1490,   max_mismatch)
leray <- find_primer(seqs, primers$mlCOIintF, max_mismatch)
# El reverso se busca DESPUÉS del primer Leray (que está dentro del barcode)
hco   <- find_primer(seqs, primers$HCO2198rc, max_mismatch, from = leray$end + 1L)

primer_table <- tibble(id = names(seqs), length = width(seqs),
                       lco_end = lco$end, leray_end = leray$end, hco_start = hco$start)
primer_table
colSums(!is.na(primer_table[, c("lco_end", "leray_end", "hco_start")]))  # cuántas tienen cada primer


# -----------------------------------------------------------------------------
# (4) Recortar las 3 regiones
# -----------------------------------------------------------------------------
cut_region <- function(seqs, from, to) {
  ok  <- !is.na(from) & !is.na(to) & to >= from
  out <- seqs[ok]
  subseq(out, start = from[ok], end = to[ok])
}

regions <- list(
  complete = seqs[width(seqs) >= min_complete],
  barcode  = cut_region(seqs, lco$end + 1L,   hco$start - 1L),
  leray    = cut_region(seqs, leray$end + 1L, hco$start - 1L)
)
sapply(regions, length)                      # nº de secuencias por región
lapply(regions, function(r) summary(width(r)))

# Comparación justa: solo queries presentes en las 3 regiones
if (same_queries) {
  common  <- Reduce(intersect, lapply(regions, names))
  regions <- lapply(regions, function(r) r[names(r) %in% common])
  cat("Queries con las 3 regiones:", length(common), "\n")
}


# -----------------------------------------------------------------------------
# (5) Selección de marco de lectura (misma lógica que orf_selection_coi12.R)
# -----------------------------------------------------------------------------
gcode <- getGeneticCode(as.character(genetic_table))

translate_frame <- function(s, frame) {
  n_cod <- pmax((width(s) - frame) %/% 3, 0)
  suppressWarnings(translate(subseq(s, start = frame + 1, width = n_cod * 3),
                             genetic.code = gcode, if.fuzzy.codon = "X",
                             no.init.codon = TRUE))
}

longest_orf <- function(s, frame) {
  p <- as.character(translate_frame(DNAStringSet(s), frame - 1))
  max(nchar(strsplit(p, "\\*")[[1]]), 0)
}

select_orf <- function(s) {
  # 5a. Stops internos en cada marco (se ignora el último codón)
  stops <- sapply(0:2, function(f) {
    p <- translate_frame(s, f)
    vcountPattern("*", subseq(p, 1, width = pmax(width(p) - 1, 0)))
  })
  stops  <- matrix(stops, ncol = 3)
  n_open <- rowSums(stops == 0)
  # 5b. Marco mayoritario entre las secuencias con un único marco limpio
  uniq     <- max.col(stops == 0)[n_open == 1]
  majority <- if (length(uniq)) as.integer(names(which.max(table(uniq)))) else 1L
  # 5c. Decisión
  frame  <- integer(length(s)); reason <- character(length(s))
  for (i in seq_along(s)) {
    open <- which(stops[i, ] == 0)
    if (length(open) == 1) {
      frame[i] <- open; reason[i] <- "unico_sin_stops"
    } else if (length(open) > 1 && majority %in% open) {
      frame[i] <- majority; reason[i] <- "empate_marco_mayoritario"
    } else if (length(open) > 1) {
      l <- sapply(open, function(f) longest_orf(s[[i]], f))
      frame[i] <- open[which.max(l)]; reason[i] <- "empate_orf_mas_largo"
    } else {
      frame[i] <- which.min(stops[i, ]); reason[i] <- "ambiguo_con_stops"
    }
  }
  tibble(id = names(s), length = width(s), frame = frame, reason = reason,
         stops_f1 = stops[, 1], stops_f2 = stops[, 2], stops_f3 = stops[, 3],
         kept = reason != "ambiguo_con_stops")
}

# Quita la 3ª base de cada codón desde el marco elegido
make_coi12 <- function(s, frame) {
  x <- substring(as.character(s), frame)
  out <- DNAStringSet(gsub("(..).", "\\1", x, perl = TRUE))
  names(out) <- names(s)
  out
}


# -----------------------------------------------------------------------------
# (6) Generar 123 y 12 de cada región y guardar
# -----------------------------------------------------------------------------
report <- list()
for (r in names(regions)) {
  s   <- regions[[r]]
  orf <- select_orf(s)
  keep <- orf$kept
  writeXStringSet(s[keep], file.path(out_dir, paste0(r, "_123.fasta")), width = 20000)
  writeXStringSet(make_coi12(s[keep], orf$frame[keep]),
                  file.path(out_dir, paste0(r, "_12.fasta")), width = 20000)
  report[[r]] <- mutate(orf, region = r, .before = 1)
  cat(sprintf("%-9s %4d secuencias | conservadas %4d | marco +%s mayoritario\n",
              r, length(s), sum(keep), names(which.max(table(orf$frame)))))
}
report <- bind_rows(report)
write_tsv(report, file.path(out_dir, "regions_report.tsv"))


# -----------------------------------------------------------------------------
# PASO 7. Revisar
# -----------------------------------------------------------------------------
report |> count(region, frame)          # complete -> +1 (empieza en el ATG); leray -> +2
report |> count(region, reason)
report |> group_by(region) |> summarise(media_pb = mean(length), min = min(length), max = max(length))

# Siguiente paso: BLAST de cada <region>_123 contra MIDORI_123 (sin recortar)
# y de cada <region>_12 contra MIDORI_12, con el mismo bloqueo taxonómico.
# Después: compare_regions_by_rank.R