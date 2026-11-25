# orf_selection_coi12.R

#  Selección del marco de lectura (ORF) de cada secuencia y generación de coi12
#  -> Adaptado a la estructura del proyecto Evaluating_coi12

#   1. Abre Evaluating_coi12.Rproj en RStudio (así here() encuentra las rutas).
#   2. En el PASO 1 eliges qué archivo procesar (Calanoida o MIDORI).

#  ENTRADAS (carpeta data/; en el PASO 1 solo cambias el nombre del archivo)
#   - all_together_calanoida_123.fasta -> en realidad es un CSV:
#       taxid,secuencia,accession,nombre_cientifico
#   - midori_123.fasta                 -> FASTA normal (>accession / secuencia)
#   El script detecta solo el formato.

#  SALIDAS (carpeta data/orf/)
#   - <nombre>_12_orf.fasta   secuencias coi12 con el marco de cada una
#   - <nombre>_orf_report.tsv informe: marco, stops por marco, motivo, kept
#   - <nombre>_distintas_vs_12.tsv  secuencias que cambian respecto a Generate_coi12.sh
#   - figuras/<nombre>/         gráficos y esquemas (PASO 10)

#  CRITERIOS BIOLÓGICOS
#   1. Código mitocondrial de invertebrados (tabla 5).
#   2. Fragmento interno codificante -> el marco correcto no tiene stops internos.
#   3. Empate (varios marcos sin stops) -> marco mayoritario del dataset;
#      si no está entre los candidatos, el ORF más largo.
#   4. Stops en los 3 marcos -> posible NUMT/error -> kept = FALSE.



# -----------------------------------------------------------------------------
# (0) Librerías
# -----------------------------------------------------------------------------
# Si falta alguna:
#   install.packages(c("here", "readr", "dplyr"))
#   install.packages("BiocManager"); BiocManager::install("Biostrings")

  library(Biostrings)
  library(here)
  library(readr)
  library(dplyr)



# -----------------------------------------------------------------------------
# (1) PARÁMETROS  <-- lo único que tienes que cambiar
# -----------------------------------------------------------------------------
# Archivos disponibles en data/ (para ver los nombres exactos):
list.files(here("data"))
#   all_together_calanoida_123.fasta  -> queries de Calanoida (CSV)
#   midori_123.fasta                  -> base de referencia MIDORI (FASTA, ~3 min)

archivo        <- "all_together_calanoida_123.fasta"   # <-- cambia aquí
# archivo      <- "midori_123.fasta"

input_file     <- here("data", archivo)
out_dir        <- here("results", "orf")
genetic_table  <- 5       # 5 = mitocondrial invertebrados (2 = vertebrados)
min_length     <- 9       # secuencias más cortas no se pueden evaluar
keep_ambiguous <- FALSE   # TRUE = conservar las que tienen stops en los 3 marcos

stopifnot("No encuentro el archivo en data/" = file.exists(input_file))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
base_name <- sub("_?123\\.fasta$", "", basename(input_file))   # p.ej. "all_together_calanoida"


# -----------------------------------------------------------------------------
# (2) Leer el archivo (detecta CSV o FASTA)
# -----------------------------------------------------------------------------
first_char <- substr(readLines(input_file, n = 1), 1, 1)

if (first_char == ">") {
  # FASTA normal (MIDORI)
  seqs <- readDNAStringSet(input_file)
} else {
  # CSV de Calanoida: taxid, seq, accession, scientific_name
  tab <- read_csv(input_file,
                  col_names = c("taxid", "seq", "accession", "scientific_name"),
                  show_col_types = FALSE)
  seqs <- DNAStringSet(tab$seq)
  # Header igual que en Preparar_inputs.Rmd: accession;nombre
  names(seqs) <- paste0(tab$accession, ";", tab$scientific_name)
}

# Mayúsculas y sin gaps; las bases raras (no IUPAC) se cambian por N
seqs <- DNAStringSet(gsub("-", "", toupper(as.character(seqs))))
cat("Secuencias leídas:", length(seqs), "\n")
summary(width(seqs))


# -----------------------------------------------------------------------------
# (3) Traducir los 3 marcos y contar stops internos (vectorizado)
# -----------------------------------------------------------------------------
# Se trabaja con todas las secuencias a la vez -> rápido incluso con MIDORI.
gcode <- getGeneticCode(as.character(genetic_table))
evaluable <- width(seqs) >= min_length

translate_frame <- function(seqs, frame) {
  w     <- width(seqs)
  n_cod <- pmax((w - frame) %/% 3, 0)                 # codones completos
  sub   <- subseq(seqs, start = frame + 1, width = n_cod * 3)
  suppressWarnings(translate(sub, genetic.code = gcode,
                             if.fuzzy.codon = "X", no.init.codon = TRUE))
}

stops <- sapply(0:2, function(f) {
  prot <- translate_frame(seqs[evaluable], f)
  # Stops INTERNOS: quitamos el último aminoácido (stop final = efecto del recorte)
  internal <- subseq(prot, start = 1, width = pmax(width(prot) - 1, 0))
  vcountPattern("*", internal)
})
stops <- matrix(stops, ncol = 3, dimnames = list(NULL, c("f1", "f2", "f3")))
head(stops)


# -----------------------------------------------------------------------------
# (4) Marco mayoritario (evidencia: secuencias con un ÚNICO marco limpio)
# -----------------------------------------------------------------------------
n_open   <- rowSums(stops == 0)
unique_f <- apply(stops == 0, 1, which)[n_open == 1]
majority <- as.integer(names(which.max(table(unlist(unique_f)))))   # 1, 2 o 3
cat("Marco mayoritario: +", majority, "\n", sep = "")
table(unlist(unique_f))


# -----------------------------------------------------------------------------
# (5) Elegir el marco de cada secuencia
# -----------------------------------------------------------------------------
# ORF más largo de un marco (solo se calcula cuando hace falta desempatar)
longest_orf <- function(s, frame) {
  p <- as.character(translate_frame(DNAStringSet(s), frame - 1))
  max(nchar(strsplit(p, "\\*")[[1]]), 0)
}

choose_frame <- function(i) {
  st   <- stops[i, ]
  open <- which(st == 0)
  if (length(open) == 1) return(c(open, 1))                     # único sin stops
  if (length(open) > 1) {
    if (majority %in% open) return(c(majority, 2))              # empate -> mayoritario
    s <- seqs[evaluable][[i]]
    l <- sapply(open, function(f) longest_orf(s, f))
    return(c(open[which.max(l)], 3))                            # empate -> ORF más largo
  }
  c(which.min(st), 4)                                           # ambiguo (stops en todos)
}

reasons <- c("unico_sin_stops", "empate_marco_mayoritario",
             "empate_orf_mas_largo", "ambiguo_con_stops")

# Rápido: la mayoría de casos se resuelven sin bucle
frame  <- integer(nrow(stops)); reason <- integer(nrow(stops))
easy   <- n_open == 1
frame[easy]  <- max.col(stops[easy, , drop = FALSE] == 0)
reason[easy] <- 1
tie_maj <- n_open > 1 & stops[, majority] == 0
frame[tie_maj] <- majority; reason[tie_maj] <- 2
rest <- which(!easy & !tie_maj)                    # pocos casos -> bucle
for (i in rest) { r <- choose_frame(i); frame[i] <- r[1]; reason[i] <- r[2] }


# -----------------------------------------------------------------------------
# (6) Informe
# -----------------------------------------------------------------------------
report <- tibble(
  id     = names(seqs),
  length = width(seqs),
  frame  = NA_integer_,
  reason = "demasiado_corta",
  stops_f1 = NA_integer_, stops_f2 = NA_integer_, stops_f3 = NA_integer_
)
report$frame[evaluable]    <- frame
report$reason[evaluable]   <- reasons[reason]
report$stops_f1[evaluable] <- stops[, 1]
report$stops_f2[evaluable] <- stops[, 2]
report$stops_f3[evaluable] <- stops[, 3]
report <- report |>
  mutate(kept = reason %in% reasons[1:3] | (keep_ambiguous & reason == reasons[4]))

# Resúmenes para revisar
count(report, frame)
count(report, reason)
report |> filter(frame != majority | reason != "unico_sin_stops")   # casos a mirar


# -----------------------------------------------------------------------------
# PASO 7. Generar coi12 con el marco de cada secuencia
# -----------------------------------------------------------------------------
# Desde el marco elegido, el patrón "(..)." quita 1 de cada 3 bases = 3ª posición.
keep_idx <- which(report$kept)
s_chr    <- as.character(seqs[keep_idx])
s_chr    <- substring(s_chr, report$frame[keep_idx])        # empezar en el codón 1
coi12    <- gsub("(..).", "\\1", s_chr, perl = TRUE)

out12 <- DNAStringSet(coi12)
names(out12) <- report$id[keep_idx]


# -----------------------------------------------------------------------------
# PASO 8. Guardar
# -----------------------------------------------------------------------------
fasta_out  <- file.path(out_dir, paste0(base_name, "_12_orf.fasta"))
report_out <- file.path(out_dir, paste0(base_name, "_orf_report.tsv"))
writeXStringSet(out12, fasta_out, width = 20000)          # una línea por secuencia
write_tsv(report, report_out)

cat("\nConservadas:", length(keep_idx), "de", nrow(report),
    "\nFASTA coi12:", fasta_out, "\nInforme:   ", report_out, "\n")


# -----------------------------------------------------------------------------
# PASO 9. Comparar con el coi12 actual (hecho con Generate_coi12.sh, marco +2)
# -----------------------------------------------------------------------------
# Busca automáticamente la pareja *_12.fasta en data/ y compara por ACCESSION
# (midori_123 usa ">AB000119_37099" y midori_12 solo ">AB000119").
file_12 <- sub("_123\\.fasta$", "_12.fasta", input_file)

read_any <- function(f) {                     # FASTA o CSV del proyecto
  if (substr(readLines(f, n = 1), 1, 1) == ">") {
    x <- readDNAStringSet(f); setNames(as.character(x), names(x))
  } else {
    t <- read_csv(f, col_names = c("taxid", "seq", "accession", "scientific_name"),
                  show_col_types = FALSE)
    setNames(t$seq, t$accession)
  }
}
acc <- function(x) sub("[_; ].*$", "", x)     # accession = antes de "_" o ";"

if (file.exists(file_12)) {
  old_12  <- read_any(file_12)
  names(old_12) <- acc(names(old_12))
  new_12  <- setNames(as.character(out12), acc(names(out12)))
  comunes <- intersect(names(new_12), names(old_12))
  iguales <- old_12[comunes] == new_12[comunes]
  cat("\nComparación con", basename(file_12), ":\n",
      " secuencias en común:", length(comunes), "\n",
      " idénticas:          ", sum(iguales), sprintf("(%.2f %%)", 100 * mean(iguales)), "\n",
      " distintas:          ", sum(!iguales), "\n")
  # Las distintas, con su marco y motivo:
  distintas <- report |> filter(acc(id) %in% comunes[!iguales])
  distintas
  write_tsv(distintas, file.path(out_dir, paste0(base_name, "_distintas_vs_12.tsv")))
} else {
  cat("No hay", basename(file_12), "en data/: se omite la comparación.\n")
}

# =============================================================================
# PASO 10. FIGURAS Y ESQUEMAS
# =============================================================================
# Se guardan en <out_dir>/figuras como PNG (300 dpi, listos para el TFG).
# Paquetes extra: install.packages(c("ggplot2", "tidyr", "patchwork"))
library(ggplot2)
library(tidyr)

fig_dir <- file.path(out_dir, "figuras", base_name)   # una carpeta por archivo
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tema <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

col_frame  <- c("+1" = "#E69F00", "+2" = "#0072B2", "+3" = "#009E73", "NA" = "grey60")
col_reason <- c(unico_sin_stops          = "#0072B2",
                empate_marco_mayoritario = "#56B4E9",
                empate_orf_mas_largo     = "#E69F00",
                ambiguo_con_stops        = "#D55E00",
                demasiado_corta          = "grey60")

report_plot <- report |>
  mutate(frame_lab = ifelse(is.na(frame), "NA", paste0("+", frame)),
         reason    = factor(reason, names(col_reason)))


# -----------------------------------------------------------------------------
# 10a. ESQUEMA: cómo se genera coi12 a partir del marco elegido
# -----------------------------------------------------------------------------
# Toma una secuencia real del dataset y dibuja sus primeras 30 bases:
# cada base coloreada por su posición en el codón; las 3ª posiciones se tachan.
ej     <- which(report$kept)[1]
ej_seq <- strsplit(substr(as.character(seqs[ej]), 1, 30), "")[[1]]
ej_fr  <- report$frame[ej]

esquema <- tibble(pos = seq_along(ej_seq), base = ej_seq) |>
  mutate(cod_pos = ifelse(pos < ej_fr, "fuera de marco",
                          paste0(((pos - ej_fr) %% 3) + 1, "ª posición")),
         codon   = ifelse(pos < ej_fr, NA, (pos - ej_fr) %/% 3 + 1),
         quitada = cod_pos %in% c("3ª posición", "fuera de marco"))

p_esquema <- ggplot(esquema, aes(pos, 1)) +
  geom_tile(aes(fill = cod_pos), colour = "white", height = 0.8) +
  geom_text(aes(label = base), fontface = "bold", size = 4.5) +
  geom_segment(data = filter(esquema, quitada),
               aes(x = pos - 0.4, xend = pos + 0.4, y = 0.65, yend = 1.35),
               colour = "red", linewidth = 0.9) +
  # llaves de codón
  geom_segment(data = esquema |> filter(!is.na(codon)) |> group_by(codon) |>
                 summarise(x = min(pos), xend = max(pos)),
               aes(x = x - 0.4, xend = xend + 0.4, y = 1.5, yend = 1.5),
               linewidth = 0.6, inherit.aes = FALSE) +
  annotate("text", x = 0.5, y = 0.3, hjust = 0, size = 3.5,
           label = paste0("Secuencia: ", report$id[ej], "  |  marco elegido: +", ej_fr,
                          "  ->  coi12: ", substr(coi12[1], 1, 20), "...")) +
  scale_fill_manual(values = c("1ª posición" = "#0072B2", "2ª posición" = "#56B4E9",
                               "3ª posición" = "#F0E442", "fuera de marco" = "grey80")) +
  scale_y_continuous(limits = c(0.2, 1.7)) +
  labs(title = "Generación de coi12: se elimina la 3ª posición de cada codón",
       subtitle = "Las líneas negras marcan los codones; las bases tachadas no pasan a coi12",
       fill = NULL, x = "Posición en la secuencia", y = NULL) +
  tema + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
               panel.grid = element_blank(), legend.position = "bottom")
p_esquema
ggsave(file.path(fig_dir, "01_esquema_coi12.png"), p_esquema, width = 11, height = 3.8, dpi = 300)


# -----------------------------------------------------------------------------
# 10b. ESQUEMA: árbol de decisión con el nº de secuencias en cada rama
# -----------------------------------------------------------------------------
n_r <- table(factor(report$reason, names(col_reason)))
cajas <- tribble(
  ~x,  ~y,  ~lab,                                                         ~tipo,
  3,   5,   paste0("Secuencias leídas\nn = ", nrow(report)),              "inicio",
  1,   4,   paste0("< ", min_length, " pb\nn = ", n_r["demasiado_corta"]), "fuera",
  3.5, 4,   "Traducir 3 marcos (tabla 5)\ny contar stops internos",       "paso",
  1.2, 2.6, paste0("1 marco sin stops\nn = ", n_r["unico_sin_stops"]),    "ok",
  3.5, 2.6, "Varios marcos\nsin stops",                                   "paso",
  5.8, 2.6, paste0("Stops en los 3 marcos\n(¿NUMT?)\nn = ", n_r["ambiguo_con_stops"]), "fuera",
  2.6, 1.2, paste0("Marco mayoritario (+", majority, ")\nn = ",
                   n_r["empate_marco_mayoritario"]),                      "ok",
  4.6, 1.2, paste0("ORF más largo\nn = ", n_r["empate_orf_mas_largo"]),   "ok"
)
flechas <- tribble(
  ~x, ~y, ~xend, ~yend,
  3, 4.65, 1, 4.35,    3, 4.65, 3.5, 4.35,
  3.5, 3.65, 1.2, 2.95, 3.5, 3.65, 3.5, 2.95, 3.5, 3.65, 5.8, 3.05,
  3.5, 2.25, 2.6, 1.55, 3.5, 2.25, 4.6, 1.55
)
p_arbol <- ggplot() +
  geom_segment(data = flechas, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(2.5, "mm")), colour = "grey40") +
  geom_label(data = cajas, aes(x, y, label = lab, fill = tipo),
             size = 3.6, label.padding = unit(3, "mm"), lineheight = 0.95,
             key_glyph = draw_key_rect) +
  scale_fill_manual(values = c(inicio = "grey90", paso = "#DCEBF7",
                               ok = "#CDEBD3", fuera = "#F6D0C4"),
                    labels = c(inicio = "Entrada", paso = "Paso", ok = "Se conserva",
                               fuera = "Se descarta")) +
  coord_cartesian(xlim = c(0.2, 6.6), ylim = c(0.6, 5.5)) +
  labs(title = "Criterio de selección del marco de lectura", fill = NULL) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
p_arbol
ggsave(file.path(fig_dir, "02_arbol_decision.png"), p_arbol, width = 10, height = 6.5, dpi = 300)


# -----------------------------------------------------------------------------
# 10c. BARRAS: marco de lectura elegido
# -----------------------------------------------------------------------------
p_frames <- report_plot |>
  count(frame_lab) |>
  mutate(pct = 100 * n / sum(n)) |>
  ggplot(aes(frame_lab, n, fill = frame_lab)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)), vjust = -0.2, size = 3.8) +
  scale_fill_manual(values = col_frame) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Marco de lectura elegido", x = "Marco", y = "Nº de secuencias") +
  tema
p_frames
ggsave(file.path(fig_dir, "03_barras_marco.png"), p_frames, width = 5.5, height = 4.5, dpi = 300)


# -----------------------------------------------------------------------------
# 10d. BARRAS: motivo de la decisión
# -----------------------------------------------------------------------------
p_reason <- report_plot |>
  count(reason, .drop = FALSE) |>
  mutate(pct = 100 * n / sum(n)) |>
  ggplot(aes(reorder(reason, n), n, fill = reason)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", n, pct)), hjust = -0.1, size = 3.6) +
  coord_flip() +
  scale_fill_manual(values = col_reason) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(title = "Motivo de la decisión", x = NULL, y = "Nº de secuencias") +
  tema
p_reason
ggsave(file.path(fig_dir, "04_barras_motivo.png"), p_reason, width = 7, height = 4, dpi = 300)


# -----------------------------------------------------------------------------
# 10e. GRÁFICO: nº de stops internos en cada marco
# -----------------------------------------------------------------------------
# En el marco correcto los stops deberían ser 0; en los otros dos, muchos.
stops_long <- report_plot |>
  filter(!is.na(stops_f1)) |>
  pivot_longer(c(stops_f1, stops_f2, stops_f3), names_to = "marco", values_to = "stops") |>
  mutate(marco = paste0("+", sub("stops_f", "", marco)))

p_stops <- ggplot(stops_long, aes(marco, stops, fill = marco)) +
  geom_violin(alpha = 0.35, colour = NA, scale = "width", show.legend = FALSE) +
  geom_boxplot(width = 0.18, outlier.size = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = col_frame) +
  labs(title = "Codones stop internos en cada marco de lectura",
       subtitle = "El marco correcto concentra las secuencias en 0 stops",
       x = "Marco", y = "Nº de stops internos") +
  tema
p_stops
ggsave(file.path(fig_dir, "05_stops_por_marco.png"), p_stops, width = 6, height = 4.5, dpi = 300)


# -----------------------------------------------------------------------------
# 10f. MAPA DE CALOR: combinaciones de marcos sin stops
# -----------------------------------------------------------------------------
# Qué marcos quedan "abiertos" (0 stops) en cada secuencia -> explica los empates.
p_patron <- report_plot |>
  filter(!is.na(stops_f1)) |>
  mutate(patron = paste0(ifelse(stops_f1 == 0, "+1 ", ""),
                         ifelse(stops_f2 == 0, "+2 ", ""),
                         ifelse(stops_f3 == 0, "+3", "")),
         patron = ifelse(trimws(patron) == "", "ninguno", trimws(patron))) |>
  count(patron, frame_lab) |>
  ggplot(aes(frame_lab, patron, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n), size = 4) +
  scale_fill_gradient(low = "#DCEBF7", high = "#08306B", trans = "log10") +
  labs(title = "Marcos sin stops vs marco elegido",
       x = "Marco elegido", y = "Marcos sin stops internos", fill = "Nº sec.") +
  tema
p_patron
ggsave(file.path(fig_dir, "06_heatmap_marcos_abiertos.png"), p_patron, width = 6, height = 4.5, dpi = 300)


# -----------------------------------------------------------------------------
# 10g. HISTOGRAMA: longitudes 123 vs coi12
# -----------------------------------------------------------------------------
long_df <- bind_rows(
  tibble(version = "123 (original)", length = width(seqs[keep_idx])),
  tibble(version = "12 (coi12)",     length = width(out12))
) |> mutate(version = factor(version, c("123 (original)", "12 (coi12)")))
p_len <- ggplot(long_df, aes(length, fill = version)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  geom_vline(data = long_df |> group_by(version) |> summarise(m = median(length)),
             aes(xintercept = m, colour = version), linetype = "dashed",
             show.legend = FALSE) +
  scale_fill_manual(values = c("123 (original)" = "#0072B2", "12 (coi12)" = "#E69F00")) +
  scale_colour_manual(values = c("123 (original)" = "#0072B2", "12 (coi12)" = "#E69F00")) +
  labs(title = "Longitud de las secuencias antes y después de quitar la 3ª base",
       subtitle = "Líneas discontinuas: mediana (coi12 ≈ 2/3 de la original)",
       x = "Longitud (pb)", y = "Nº de secuencias", fill = NULL) +
  tema + theme(legend.position = "top")
p_len
ggsave(file.path(fig_dir, "07_histograma_longitudes.png"), p_len, width = 8, height = 4.5, dpi = 300)


# -----------------------------------------------------------------------------
# 10h. PANEL RESUMEN (todas juntas, opcional con patchwork)
# -----------------------------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  panel <- (p_frames | p_reason) / (p_stops | p_patron) / p_len +
    plot_annotation(title = paste("Selección de ORF:", base_name),
                    tag_levels = "A",
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
  ggsave(file.path(fig_dir, "00_panel_resumen.png"), panel, width = 13, height = 13, dpi = 300)
}

cat("Figuras guardadas en:", fig_dir, "\n")
list.files(fig_dir)