# =============================================================================
#  compare_regions_by_rank.R
#  Compara la precisión de asignación taxonómica de las 3 regiones de COI
#  (complete, barcode, leray) en sus versiones 123 y 12, según el rango
#  taxonómico bloqueado en el BLAST.
# =============================================================================
#
#  ENTRADA (carpeta in_dir): un CSV por región y versión, con el formato que ya
#  genera analyze_blast_results.R (individual_results, tras el LCA):
#       join_complete_123.csv   join_complete_12.csv
#       join_barcode_123.csv    join_barcode_12.csv
#       join_leray_123.csv      join_leray_12.csv
#  Columnas que se usan: qseqid, blocked_rank, pident, precision
#  (si un archivo tiene varias filas por query y rango -los top-5 hits- se queda
#   con el mejor hit por pident y avisa).
#
#  QUÉ CALCULA
#   1. % de queries bien asignadas HASTA cada nivel (kingdom ... species).
#   2. % de queries asignadas al nivel MÁXIMO POSIBLE dado el bloqueo:
#        no_block -> species | species -> genus | genus -> family
#        family   -> order   | order   -> class | class -> phylum
#   3. Test de McNemar exacto (pareado) 123 vs 12 dentro de cada región y rango:
#      ¿quitar la 3ª base cambia el acierto al nivel máximo posible?
#
#  SALIDAS (out_dir): resumen_por_nivel.csv, resumen_max_posible.csv,
#                     mcnemar_123_vs_12.csv, fig_por_nivel.png, fig_max_posible.png
# =============================================================================


# -----------------------------------------------------------------------------
# PASO 0. Librerías
# -----------------------------------------------------------------------------
suppressMessages({
  library(tidyverse)
  library(here)
})


# -----------------------------------------------------------------------------
# PASO 1. PARÁMETROS
# -----------------------------------------------------------------------------
in_dir  <- here("05_blast_results", "regions")
out_dir <- here("08_orf_results", "comparacion_regiones")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

levels_tax  <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
blocked_ord <- c("no_block", "species", "genus", "family", "order", "class")
max_possible <- c(no_block = "species", species = "genus", genus = "family",
                  family = "order", order = "class", class = "phylum")
region_ord  <- c("complete", "barcode", "leray")


# -----------------------------------------------------------------------------
# PASO 2. Leer los 6 archivos
# -----------------------------------------------------------------------------
files <- list.files(in_dir, pattern = "^join_(complete|barcode|leray)_(123|12)\\.csv$",
                    full.names = TRUE)
stopifnot("No hay archivos join_<region>_<version>.csv en in_dir" = length(files) > 0)

read_one <- function(f) {
  parts <- str_match(basename(f), "^join_(\\w+)_(123|12)\\.csv$")
  read_csv(f, col_types = cols(.default = col_character()), show_col_types = FALSE) |>
    filter(blocked_rank != "blocked_rank") |>          # cabeceras repetidas al concatenar
    mutate(region = parts[2], version = parts[3], pident = as.numeric(pident))
}
res <- map_dfr(files, read_one)

# Una fila por query x rango bloqueado (si hay varias, el mejor hit)
dups <- res |> count(region, version, qseqid, blocked_rank) |> filter(n > 1)
if (nrow(dups) > 0) {
  warning("Hay varias filas por query y rango (¿top-5 hits?). ",
          "Se usa el mejor hit por pident. Lo ideal es usar los resultados del LCA.")
}
res <- res |>
  arrange(desc(pident)) |>
  distinct(region, version, qseqid, blocked_rank, .keep_all = TRUE)

# Nivel alcanzado como número: 0 = nada (Worse, Unassigned, No_Hit), 1 = kingdom ... 7 = species
res <- res |>
  mutate(level_num  = match(precision, levels_tax, nomatch = 0),
         max_num    = match(max_possible[blocked_rank], levels_tax),
         at_max     = level_num >= max_num,
         region     = factor(region, region_ord),
         version    = factor(version, c("123", "12")),
         blocked_rank = factor(blocked_rank, blocked_ord))

count(res, region, version)   # nº de filas por combinación


# -----------------------------------------------------------------------------
# PASO 3. % bien asignadas hasta cada nivel
# -----------------------------------------------------------------------------
by_level <- res |>
  crossing(level = factor(levels_tax, levels_tax)) |>
  group_by(region, version, blocked_rank, level) |>
  summarise(n = n(),
            pct_correct = 100 * mean(level_num >= as.integer(level)),
            .groups = "drop")
write_csv(by_level, file.path(out_dir, "resumen_por_nivel.csv"))


# -----------------------------------------------------------------------------
# PASO 4. % asignadas al nivel máximo posible
# -----------------------------------------------------------------------------
by_max <- res |>
  group_by(region, version, blocked_rank) |>
  summarise(n = n(),
            nivel_max = max_possible[as.character(first(blocked_rank))],
            pct_max_posible = 100 * mean(at_max),
            .groups = "drop")
write_csv(by_max, file.path(out_dir, "resumen_max_posible.csv"))
by_max |>
  select(region, version, blocked_rank, pct_max_posible) |>
  pivot_wider(names_from = blocked_rank, values_from = pct_max_posible)


# -----------------------------------------------------------------------------
# PASO 5. McNemar exacto pareado: 123 vs 12 (mismas queries)
# -----------------------------------------------------------------------------
# b = aciertos solo en 123 ; c = aciertos solo en 12. Con b + c pares
# discordantes, el test exacto es un binomial(b, b + c, 0.5).
mcnemar <- res |>
  select(region, version, qseqid, blocked_rank, at_max) |>
  pivot_wider(names_from = version, values_from = at_max, names_prefix = "v") |>
  filter(!is.na(v123), !is.na(v12)) |>
  group_by(region, blocked_rank) |>
  summarise(n_pares   = n(),
            solo_123  = sum(v123 & !v12),
            solo_12   = sum(!v123 & v12),
            .groups = "drop") |>
  mutate(p_valor = map2_dbl(solo_123, solo_12, \(b, c)
           if (b + c == 0) 1 else binom.test(b, b + c, 0.5)$p.value),
         p_ajustado_BH = p.adjust(p_valor, method = "BH"))
write_csv(mcnemar, file.path(out_dir, "mcnemar_123_vs_12.csv"))
mcnemar


# -----------------------------------------------------------------------------
# PASO 6. Figuras
# -----------------------------------------------------------------------------
# 6a. Curvas: % bien asignado hasta cada nivel, por rango bloqueado
p1 <- ggplot(by_level, aes(level, pct_correct, colour = region,
                           linetype = version, group = interaction(region, version))) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
  facet_wrap(~ blocked_rank, nrow = 2) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Nivel taxonómico", y = "% queries bien asignadas hasta ese nivel",
       colour = "Región", linetype = "Versión",
       title = "Precisión de asignación por región de COI y rango bloqueado") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "fig_por_nivel.png"), p1, width = 11, height = 6.5, dpi = 300)

# 6b. Barras: % al nivel máximo posible
p2 <- ggplot(by_max, aes(blocked_rank, pct_max_posible, fill = version)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_wrap(~ region) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Rango bloqueado", y = "% asignadas al nivel máximo posible",
       fill = "Versión", title = "123 vs 12 en cada región de COI") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "fig_max_posible.png"), p2, width = 10, height = 4.5, dpi = 300)

p1; p2