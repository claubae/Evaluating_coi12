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

#  SALIDAS (carpeta 08_orf_results/<nombre>/)
#   - <nombre>_12_orf.fasta   secuencias coi12 con el marco de cada una
#   - <nombre>_orf_report.tsv informe: marco, stops por marco, motivo, kept
#   - <nombre>_distintas_vs_12.tsv  secuencias que cambian respecto a Generate_coi12.sh
#   - figuras/                  gráficos y esquemas (PASO 10) y estadística (PASO 11)
#   - estadisticas.csv          proporciones con IC 95 % (PASO 11)

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
base_name      <- sub("_?123\\.fasta$", "", archivo)       # p.ej. "all_together_calanoida"
out_dir        <- here("08_orf_results", base_name)        # una carpeta por archivo
genetic_table  <- 5       # 5 = mitocondrial invertebrados (2 = vertebrados)
min_length     <- 9       # secuencias más cortas no se pueden evaluar
keep_ambiguous <- FALSE   # TRUE = conservar las que tienen stops en los 3 marcos

stopifnot("No encuentro el archivo en data/" = file.exists(input_file))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


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

fig_dir <- file.path(out_dir, "figuras")
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
       subtitle = "Líneas discontinuas: mediana (coi12 = aprox. 2/3 de la original)",
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

# =============================================================================
# PASO 11. ESTADÍSTICA Y FIGURAS DE RESULTADOS (para presentar)
# =============================================================================
# Responde a tres preguntas:
#   1. ¿Qué proporción de secuencias NO está en el marco asumido (+2)?
#      -> proporciones con intervalo de confianza del 95 % (Wilson)
#   2. ¿Depende de la longitud de la secuencia?
#   3. ¿Cuánto cambia el coi12 de las secuencias corregidas respecto al de
#      Generate_coi12.sh? -> identidad entre coi12 antiguo y corregido
# Referencia del IC de Wilson: Brown, Cai & DasGupta (2001) Statistical Science
#   16:101-133. https://doi.org/10.1214/ss/1009213286

# -----------------------------------------------------------------------------
# 11a. Proporciones con IC 95 % (Wilson)
# -----------------------------------------------------------------------------
wilson <- function(x, n, z = 1.96) {
  p <- x / n; d <- 1 + z^2 / n
  centro <- (p + z^2 / (2 * n)) / d
  margen <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(pct = 100 * p, ic_inf = 100 * max(0, centro - margen), ic_sup = 100 * min(1, centro + margen))
}

ev_rep <- report |> filter(reason != "demasiado_corta")
n_ev   <- nrow(ev_rep)
categorias <- list(
  "Marco +2 (el asumido)"            = ev_rep$frame == 2 & ev_rep$kept,
  "Marco +1 (corregidas)"            = ev_rep$frame == 1 & ev_rep$kept,
  "Marco +3 (corregidas)"            = ev_rep$frame == 3 & ev_rep$kept,
  "Fuera de +2 (total corregidas)"   = ev_rep$frame != 2 & ev_rep$kept,
  "Excluidas: stops en los 3 marcos" = ev_rep$reason == "ambiguo_con_stops",
  "Resueltas por desempate"          = ev_rep$reason %in% c("empate_marco_mayoritario",
                                                            "empate_orf_mas_largo")
)
estadisticas <- bind_rows(lapply(names(categorias), function(k) {
  x <- sum(categorias[[k]])
  tibble(categoria = k, n = x, total = n_ev, !!!as.list(wilson(x, n_ev)),
         por_millon = 1e6 * x / n_ev)
}))
if (exists("comunes")) {
  estadisticas <- bind_rows(estadisticas, tibble(
    categoria = "Idénticas a la versión 12 actual", n = sum(iguales), total = length(comunes),
    !!!as.list(wilson(sum(iguales), length(comunes))),
    por_millon = 1e6 * sum(iguales) / length(comunes)))
}
estadisticas
write_csv(estadisticas, file.path(out_dir, "estadisticas.csv"))


# -----------------------------------------------------------------------------
# 11b. FIGURA: proporción de cada categoría con IC 95 % (escala logarítmica)
# -----------------------------------------------------------------------------
# Escala log porque conviven 99,99 % y 0,002 %: así se ven todas.
est_plot <- estadisticas |>
  filter(n > 0, !grepl("Idénticas", categoria)) |>
  mutate(categoria = factor(categoria, rev(categoria)),
         grupo = ifelse(grepl("\\+2 \\(el", categoria), "Asumido", 
                 ifelse(grepl("Excluidas", categoria), "Excluidas", "Corregidas / desempate")))
p_ic <- ggplot(est_plot, aes(pct, categoria, colour = grupo)) +
  geom_errorbar(aes(xmin = pmax(ic_inf, 1e-4), xmax = ic_sup), width = 0.25,
                linewidth = 0.8, orientation = "y") +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("%s sec. · %s %%",
                                format(n, big.mark = ".", decimal.mark = ","),
                                formatC(signif(pct, ifelse(pct > 99, 6, 2)), format = "fg", decimal.mark = ","))),
            vjust = -1.1, size = 3.4, show.legend = FALSE) +
  scale_x_log10(labels = function(x) paste0(format(x, decimal.mark = ",", scientific = FALSE,
                                                   drop0trailing = TRUE), " %"),
                expand = expansion(mult = 0.18)) +
  coord_cartesian(clip = "off") +
  scale_colour_manual(values = c("Asumido" = "#0072B2", "Corregidas / desempate" = "#E69F00",
                                 "Excluidas" = "#D55E00")) +
  labs(title = "¿Cuántas secuencias NO están en el marco +2?",
       subtitle = sprintf("%s secuencias evaluadas · punto = proporción, barra = IC 95 %% (Wilson)",
                          format(n_ev, big.mark = ".", decimal.mark = ",")),
       x = "% de secuencias (escala logarítmica)", y = NULL, colour = NULL) +
  tema + theme(legend.position = "bottom")
p_ic
ggsave(file.path(fig_dir, "08_proporciones_IC95.png"), p_ic, width = 10, height = 5.5, dpi = 300)


# -----------------------------------------------------------------------------
# 11c. FIGURA: ¿las secuencias fuera de +2 tienen longitudes raras?
# -----------------------------------------------------------------------------
# La región Leray mide ~313 pb. Si una secuencia tiene 1-2 pb de más o de menos
# al inicio (recorte desigual), su marco se desplaza.
len_bins <- c(0, 200, 300, 310, 312, 313, 314, 316, 330, Inf)
len_lab  <- c("<200", "200-299", "300-309", "310-311", "312", "313", "314-315", "316-329", ">=330")
por_long <- ev_rep |>
  filter(kept) |>
  mutate(rango = cut(length, len_bins, labels = len_lab, right = FALSE)) |>
  group_by(rango) |>
  summarise(n = n(), fuera = sum(frame != 2), .groups = "drop") |>
  filter(n >= 30) |>            # rangos con < 30 secuencias darían IC enormes
  rowwise() |>
  mutate(w = list(wilson(fuera, n))) |>
  mutate(pct = w[["pct"]], ic_inf = w[["ic_inf"]], ic_sup = w[["ic_sup"]]) |>
  ungroup() |> select(-w)
write_csv(por_long, file.path(out_dir, "fuera_de_marco_por_longitud.csv"))

p_long <- ggplot(por_long, aes(rango, pct)) +
  geom_col(fill = "#E69F00", alpha = 0.85, width = 0.7) +
  geom_errorbar(aes(ymin = ic_inf, ymax = ic_sup), width = 0.2, colour = "grey30") +
  geom_text(aes(y = ic_sup, label = sprintf("%s\nde %s", trimws(format(fuera, big.mark = ".", decimal.mark = ",")),
                                            trimws(format(n, big.mark = ".", decimal.mark = ",")))),
            vjust = -0.3, size = 2.9, lineheight = 0.9) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Secuencias fuera del marco +2 según su longitud",
       subtitle = "Etiqueta: nº fuera de +2 de nº total del rango (rangos con 30 o más sec.) · barra = IC 95 % (Wilson)",
       x = "Longitud de la secuencia (pb)", y = "% fuera del marco +2") +
  tema
p_long
ggsave(file.path(fig_dir, "09_fuera_de_marco_por_longitud.png"), p_long, width = 9, height = 5, dpi = 300)


# -----------------------------------------------------------------------------
# 11d. FIGURA: ¿cuánto se parecía el coi12 antiguo al corregido?
# -----------------------------------------------------------------------------
# Para las secuencias corregidas comparamos posición a posición el coi12 de
# Generate_coi12.sh con el nuevo. Si el marco estaba mal, el parecido es bajo
# (se quitaron bases equivocadas). Como referencia, se compara también con un
# grupo de secuencias que NO cambian (identidad = 100 %).
identidad <- function(a, b) {
  L <- min(nchar(a), nchar(b))
  if (L == 0) return(NA_real_)
  100 * mean(strsplit(substr(a, 1, L), "")[[1]] == strsplit(substr(b, 1, L), "")[[1]])
}
p_ident <- NULL
if (exists("comunes") && any(!iguales)) {
  dif <- comunes[!iguales]
  ident <- tibble(accession = dif,
                  identidad = mapply(identidad, old_12[dif], new_12[dif]),
                  marco = report$frame[match(dif, acc(report$id))]) |>
    mutate(marco = paste0("Marco real +", marco))
  write_csv(ident, file.path(out_dir, "identidad_coi12_antiguo_vs_corregido.csv"))

  p_ident <- ggplot(ident, aes(identidad, fill = marco)) +
    geom_histogram(binwidth = 2.5, boundary = 0, colour = "white", alpha = 0.9) +
    geom_vline(xintercept = median(ident$identidad), linetype = "dashed") +
    annotate("label", x = median(ident$identidad), y = Inf, vjust = 1.3, size = 3.4,
             label = sprintf("mediana = %.1f %%", median(ident$identidad))) +
    scale_fill_manual(values = c("Marco real +1" = "#E69F00", "Marco real +3" = "#009E73",
                                 "Marco real +2" = "#0072B2")) +
    scale_x_continuous(limits = c(0, 100)) +
    labs(title = "Parecido entre el coi12 de Generate_coi12.sh y el corregido",
         subtitle = sprintf("%d secuencias corregidas · las idénticas (100 %%) son %s",
                            length(dif), format(sum(iguales), big.mark = ".", decimal.mark = ",")),
         x = "Identidad posición a posición (%)", y = "Nº de secuencias", fill = NULL) +
    tema + theme(legend.position = "top")
  p_ident
  ggsave(file.path(fig_dir, "10_identidad_coi12_antiguo_vs_corregido.png"), p_ident,
         width = 8, height = 4.8, dpi = 300)
} else {
  cat("No hay secuencias corregidas: se omite la figura 10.\n")
}


# -----------------------------------------------------------------------------
# 11e. INFOGRAFÍA RESUMEN (una sola imagen para enseñar)
# -----------------------------------------------------------------------------
fmt <- function(x) format(x, big.mark = ".", decimal.mark = ",")
get_pct <- function(k) estadisticas$pct[estadisticas$categoria == k]
tarjetas <- tibble(
  x = c(1, 2, 3, 1, 2, 3), y = c(2, 2, 2, 1, 1, 1),
  numero = c(fmt(nrow(report)),
             paste0(format(round(get_pct("Marco +2 (el asumido)"), 3), decimal.mark = ","), " %"),
             fmt(sum(ev_rep$frame != 2 & ev_rep$kept)),
             fmt(sum(ev_rep$reason == "ambiguo_con_stops")),
             fmt(sum(ev_rep$reason %in% c("empate_marco_mayoritario", "empate_orf_mas_largo"))),
             if (exists("comunes")) paste0(format(round(100 * mean(iguales), 2), decimal.mark = ","), " %") else "—"),
  texto = c("secuencias analizadas",
            "en el marco +2\n(la suposición del TFG)",
            "corregidas\n(estaban en +1 o +3)",
            "excluidas\n(stops en los 3 marcos, ¿NUMT?)",
            "resueltas por desempate\n(varios marcos sin stops)",
            "idénticas al coi12 de\nGenerate_coi12.sh"),
  color = c("#4D4D4D", "#0072B2", "#E69F00", "#D55E00", "#56B4E9", "#009E73")
)
p_info <- ggplot(tarjetas) +
  geom_tile(aes(x, y, fill = color), width = 0.94, height = 0.9, alpha = 0.12) +
  geom_tile(aes(x, y, colour = color), width = 0.94, height = 0.9, fill = NA, linewidth = 1.2) +
  geom_text(aes(x, y + 0.13, label = numero, colour = color), size = 10, fontface = "bold") +
  geom_text(aes(x, y - 0.2, label = texto), size = 4, colour = "grey25", lineheight = 0.95) +
  scale_fill_identity() + scale_colour_identity() +
  labs(title = paste("Selección del marco de lectura ·", base_name),
       subtitle = "Código mitocondrial de invertebrados (tabla 5) · marco = el que no tiene codones stop internos",
       caption = "Evaluating_coi12 · orf_seleccion_coi12.R") +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
        plot.caption = element_text(colour = "grey50"),
        plot.margin = margin(15, 15, 10, 15))
p_info
ggsave(file.path(fig_dir, "11_infografia_resumen.png"), p_info, width = 12, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 11f. PANEL DE ESTADÍSTICA (todo junto)
# -----------------------------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  panel_est <- if (!is.null(p_ident)) (p_info) / (p_ic | p_long) / p_ident else p_info / (p_ic | p_long)
  panel_est <- panel_est + plot_annotation(tag_levels = "A")
  ggsave(file.path(fig_dir, "00b_panel_estadistica.png"), panel_est,
         width = 16, height = if (!is.null(p_ident)) 19 else 13, dpi = 200)
}

cat("\nResultados en:", out_dir, "\n"); print(list.files(out_dir, recursive = TRUE))