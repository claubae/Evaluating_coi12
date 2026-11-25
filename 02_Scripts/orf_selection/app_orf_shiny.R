# =============================================================================
#  app_orf_shiny.R  —  Prototipo de interfaz web (R / Shiny) para selección de
#                      ORF y generación de coi12
# =============================================================================
#
#  QUÉ ES
#  La versión en R de app_orf_gradio.py. Permite:
#     1. Elegir un archivo de la carpeta data/ del proyecto, subir uno o pegar
#        secuencias de COI
#     2. Elegir el código genético y las opciones
#     3. Ver el marco elegido para cada secuencia y por qué
#     4. Ver gráficos y un esquema codón a codón
#     5. Comparar con la versión 12 hecha con Generate_coi12.sh (si existe)
#     6. Descargar el FASTA coi12 y el informe
#
#  Usa la MISMA lógica que orf_seleccion_coi12.R y orf.py (mismos resultados).
#
#  ARCHIVOS QUE RECONOCE EN data/ (tal y como los tienes)
#     all_together_calanoida_123.fasta  -> CSV: taxid,secuencia,accession,especie
#     all_together_calanoida_12.fasta   -> su versión 12 (para comparar)
#     midori_123.fasta                  -> FASTA (>accession_taxid)
#     midori_12.fasta                   -> su versión 12 (para comparar)
#  Si eliges un archivo *_123.fasta y existe su pareja *_12.fasta, la app
#  compara automáticamente (por accession) con lo que hizo Generate_coi12.sh.
#
#  CÓMO SE EJECUTA
#     1. Abre Evaluating_coi12.Rproj en RStudio
#     2. install.packages(c("shiny", "bslib", "DT", "here"))
#        BiocManager::install("Biostrings")
#     3. Abre este archivo y pulsa "Run App" (o: shiny::runApp("02_Scripts/app_orf_shiny.R"))
#
#  DOCUMENTACIÓN
#  - Shiny, tutorial oficial:  https://shiny.posit.co/r/getstarted/
#  - Wickham (2021) Mastering Shiny. O'Reilly. Libro gratis: https://mastering-shiny.org
#  - Biostrings: https://bioconductor.org/packages/Biostrings
#  - Fundamento biológico y referencias del método: cabecera de orf.py
#
#  NOTA: con midori_123.fasta (1,5 M secuencias) funciona, pero tarda unos
#  minutos; para trabajo rutinario es mejor orf_seleccion_coi12.R u orf.py.
# =============================================================================


# -----------------------------------------------------------------------------
# PASO 0. Librerías
# -----------------------------------------------------------------------------
suppressMessages({
  library(shiny)
  library(bslib)        # tema visual
  library(DT)           # tablas interactivas
  library(here)         # rutas relativas al proyecto
  library(Biostrings)   # leer FASTA y traducir
})
options(shiny.maxRequestSize = 1000 * 1024^2)   # permite subir hasta 1 GB

data_dir <- here("data")


# -----------------------------------------------------------------------------
# PASO 1. Constantes
# -----------------------------------------------------------------------------
genetic_codes <- c(                                    # tablas del NCBI
  "5 · Mitocondrial de invertebrados (Calanoida, artrópodos)" = 5,
  "2 · Mitocondrial de vertebrados"                          = 2,
  "4 · Mitocondrial de mohos, protozoos y celentéreos"       = 4,
  "1 · Estándar (nuclear)"                                   = 1
)
col_frame  <- c("1" = "#E69F00", "2" = "#0072B2", "3" = "#009E73")
col_reason <- c(unico_sin_stops = "#0072B2", empate_marco_mayoritario = "#56B4E9",
                empate_orf_mas_largo = "#E69F00", ambiguo_con_stops = "#D55E00",
                demasiado_corta = "grey60")


# -----------------------------------------------------------------------------
# PASO 2. Leer secuencias (FASTA o el CSV del proyecto)
# -----------------------------------------------------------------------------
read_sequences <- function(path) {
  if (substr(readLines(path, n = 1), 1, 1) == ">") {
    s <- readDNAStringSet(path)
  } else {
    tab <- read.csv(path, header = FALSE, colClasses = "character",
                    col.names = c("taxid", "seq", "accession", "scientific_name"))
    s <- DNAStringSet(tab$seq)
    names(s) <- paste0(tab$accession, ";", tab$scientific_name)
  }
  DNAStringSet(gsub("-", "", toupper(as.character(s))))
}

# Accession = lo que va antes de "_" o ";" (sirve para emparejar 123 y 12)
accession <- function(x) sub("[_; ].*$", "", x)


# -----------------------------------------------------------------------------
# PASO 3. Selección de ORF (misma lógica que orf_seleccion_coi12.R)
# -----------------------------------------------------------------------------
translate_frame <- function(seqs, frame, gcode) {
  n_cod <- pmax((width(seqs) - frame) %/% 3, 0)
  sub <- subseq(seqs, start = frame + 1, width = n_cod * 3)
  suppressWarnings(translate(sub, genetic.code = gcode,
                             if.fuzzy.codon = "X", no.init.codon = TRUE))
}

select_orf <- function(seqs, table = 5, min_length = 9, keep_ambiguous = FALSE,
                       progress = NULL) {
  gcode <- getGeneticCode(as.character(table))
  ev    <- width(seqs) >= min_length

  # 3a. Stops internos en los 3 marcos (vectorizado)
  stops <- sapply(0:2, function(f) {
    if (!is.null(progress)) progress(f + 1)
    p <- translate_frame(seqs[ev], f, gcode)
    vcountPattern("*", subseq(p, 1, width = pmax(width(p) - 1, 0)))
  })
  stops <- matrix(stops, ncol = 3)

  # 3b. Marco mayoritario (solo secuencias con un único marco limpio)
  n_open   <- rowSums(stops == 0)
  uniq     <- max.col(stops == 0, ties.method = "first")[n_open == 1]
  majority <- if (length(uniq)) as.integer(names(which.max(table(uniq)))) else 2L

  # 3c. Decisión
  frame  <- integer(nrow(stops)); reason <- character(nrow(stops))
  easy   <- n_open == 1
  frame[easy]  <- max.col(stops[easy, , drop = FALSE] == 0); reason[easy] <- "unico_sin_stops"
  tie    <- n_open > 1 & stops[, majority] == 0
  frame[tie]   <- majority; reason[tie] <- "empate_marco_mayoritario"
  for (i in which(!easy & !tie)) {
    open <- which(stops[i, ] == 0)
    if (length(open) > 1) {                         # empate sin el mayoritario
      s <- seqs[ev][i]
      l <- sapply(open, function(f) max(nchar(strsplit(
             as.character(translate_frame(s, f - 1, gcode)), "\\*")[[1]]), 0))
      frame[i] <- open[which.max(l)]; reason[i] <- "empate_orf_mas_largo"
    } else {                                        # stops en los 3 marcos
      frame[i] <- which.min(stops[i, ]); reason[i] <- "ambiguo_con_stops"
    }
  }

  # 3d. Informe
  rep <- data.frame(id = names(seqs), length = width(seqs), frame = NA_integer_,
                    reason = "demasiado_corta", stops_f1 = NA_integer_,
                    stops_f2 = NA_integer_, stops_f3 = NA_integer_)
  rep$frame[ev] <- frame; rep$reason[ev] <- reason
  rep$stops_f1[ev] <- stops[, 1]; rep$stops_f2[ev] <- stops[, 2]; rep$stops_f3[ev] <- stops[, 3]
  rep$kept <- rep$reason %in% c("unico_sin_stops", "empate_marco_mayoritario",
                                "empate_orf_mas_largo") |
              (keep_ambiguous & rep$reason == "ambiguo_con_stops")

  # 3e. coi12: desde el marco elegido, "(..)." quita la 3ª base de cada codón
  k     <- which(rep$kept)
  coi12 <- gsub("(..).", "\\1", substring(as.character(seqs[k]), rep$frame[k]), perl = TRUE)
  names(coi12) <- rep$id[k]

  list(report = rep, coi12 = coi12, majority = majority)
}


# -----------------------------------------------------------------------------
# PASO 4. Esquema codón a codón (HTML)
# -----------------------------------------------------------------------------
codon_html <- function(name, seq, frame, table, n = 90) {
  gcode <- getGeneticCode(as.character(table))
  b <- strsplit(substr(seq, 1, n), "")[[1]]
  pos_col <- c("#0072B2", "#56B4E9", "#F0E442")
  cells <- vapply(seq_along(b), function(i) {
    if (i < frame) {
      sty <- "background:#ddd;color:#555;text-decoration:line-through"
      gap <- ""
    } else {
      k   <- (i - frame) %% 3
      sty <- sprintf("background:%s;color:%s;text-decoration:%s", pos_col[k + 1],
                     ifelse(k == 2, "#000", "#fff"), ifelse(k == 2, "line-through", "none"))
      gap <- ifelse(k == 0, "margin-left:4px;", "")
    }
    sprintf("<wbr><span style='%s;%stext-decoration-color:red;text-decoration-thickness:2px;padding:2px 3px;font-family:monospace;font-weight:bold;border-radius:2px'>%s</span>",
            sty, gap, b[i])
  }, "")
  prots <- vapply(0:2, function(f) {
    p <- substr(as.character(translate_frame(DNAStringSet(seq), f, gcode)), 1, n %/% 3)
    p <- gsub("\\*", "<span style='color:red;font-weight:bold'>*</span>", p)
    sprintf("<div style='font-family:monospace'><b>+%d</b>&nbsp; %s%s</div>", f + 1, p,
            ifelse(f + 1 == frame, " <b style='color:#0072B2'>← elegido</b>", ""))
  }, "")
  coi12 <- gsub("(..).", "\\1", substring(seq, frame), perl = TRUE)
  HTML(sprintf(
    "<div style='line-height:2.1;word-break:break-all'><b>%s</b> — marco elegido <b>+%d</b> (primeras %d pb de %d)<br>%s<br><small>
     <span style='background:#0072B2;color:#fff;padding:1px 5px'>1ª</span>
     <span style='background:#56B4E9;color:#fff;padding:1px 5px'>2ª</span>
     <span style='background:#F0E442;padding:1px 5px'>3ª (se elimina)</span>
     <span style='background:#ddd;padding:1px 5px'>fuera de marco</span></small><br><br>
     <b>Traducción de los 3 marcos</b>%s<br><b>coi12:</b> <code>%s…</code></div>",
    name, frame, min(n, nchar(seq)), nchar(seq), paste(cells, collapse = ""),
    paste(prots, collapse = ""), substr(coi12, 1, 60)))
}


# -----------------------------------------------------------------------------
# PASO 5. Interfaz (UI)
# -----------------------------------------------------------------------------
data_files <- if (dir.exists(data_dir)) list.files(data_dir, pattern = "\\.(fasta|fa|csv)$") else character(0)

ui <- page_sidebar(
  title = "🧬 Selección del marco de lectura y coi12",
  theme = bs_theme(bootswatch = "flatly"),

  # ---- Barra lateral: entradas ----
  sidebar = sidebar(width = 360,
    radioButtons("source", "Origen de las secuencias",
                 c("Carpeta data/ del proyecto" = "data", "Subir archivo" = "upload",
                   "Pegar secuencias" = "paste")),
    conditionalPanel("input.source == 'data'",
      selectInput("data_file", "Archivo en data/",
                  choices = data_files,
                  selected = grep("calanoida_123", data_files, value = TRUE)[1]),
      helpText("Los *_123 se comparan automáticamente con su *_12.")),
    conditionalPanel("input.source == 'upload'",
      fileInput("upload", "FASTA o CSV (taxid,seq,accession,especie)")),
    conditionalPanel("input.source == 'paste'",
      textAreaInput("pasted", "Secuencias en FASTA", rows = 6, placeholder = ">id\nACGT...")),
    selectInput("code", "Código genético", genetic_codes),
    sliderInput("min_len", "Longitud mínima (pb) para evaluar", 3, 300, 9),
    checkboxInput("keep_amb", "Conservar secuencias con stops en los 3 marcos"),
    actionButton("run", "▶ Analizar", class = "btn-primary"),
    hr(),
    downloadButton("dl_fasta", "Descargar coi12 (FASTA)"),
    downloadButton("dl_report", "Descargar informe (TSV)")
  ),

  # ---- Panel principal: resultados ----
  navset_card_tab(
    nav_panel("Resumen",
      layout_columns(col_widths = c(5, 7),
        card(card_header("Resultado"), uiOutput("summary")),
        card(card_header("Marco elegido"), plotOutput("p_frames", height = 280))),
      layout_columns(col_widths = c(6, 6),
        card(card_header("Motivo de la decisión"), plotOutput("p_reasons", height = 280)),
        card(card_header("Stops internos por marco"), plotOutput("p_stops", height = 280)))),
    nav_panel("Esquema codón a codón",
      selectInput("seq_pick", "Secuencia (primero las que NO están en el marco mayoritario)",
                  choices = NULL, width = "100%"),
      uiOutput("scheme")),
    nav_panel("Informe", DTOutput("table")),
    nav_panel("Cómo funciona", markdown("
## ¿Cómo decide el marco?
1. **Traduce los 3 marcos** con el código genético elegido (por defecto tabla 5;
   [códigos del NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).
2. **Cuenta los codones stop internos.** En el marco correcto, un fragmento interno de COI no tiene stops.
3. **Decide:** un solo marco limpio → ese; varios → el **mayoritario** del conjunto (o el ORF más largo);
   stops en los 3 → se excluye (posible **NUMT**, [Song et al. 2008](https://doi.org/10.1073/pnas.0803076105)).
4. **Genera coi12** quitando la 3ª posición de cada codón desde ese marco.

## Referencias
- Leray et al. (2013) *Front. Zool.* 10:34. [doi](https://doi.org/10.1186/1742-9994-10-34)
- Leray et al. (2022) *Environmental DNA* 4:894–907 (MIDORI2). [doi](https://doi.org/10.1002/edn3.303)
- Pagès et al. Biostrings (Bioconductor). [web](https://bioconductor.org/packages/Biostrings)
- Nugent et al. (2020) *Genome* 63:291, paquete *coil*. [GitHub](https://github.com/CNuge/coil)
- Wickham (2021) *Mastering Shiny*. [libro gratis](https://mastering-shiny.org)
"))
  )
)


# -----------------------------------------------------------------------------
# PASO 6. Lógica del servidor
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  # 6a. Análisis: se ejecuta al pulsar "Analizar"
  res <- eventReactive(input$run, {
    path <- switch(input$source,
      data   = file.path(data_dir, input$data_file),
      upload = { req(input$upload); input$upload$datapath },
      paste  = { req(nzchar(input$pasted)); f <- tempfile(fileext = ".fasta")
                 txt <- input$pasted
                 if (!startsWith(trimws(txt), ">")) txt <- paste0(">secuencia_1\n", txt)
                 writeLines(txt, f); f })

    withProgress(message = "Analizando…", value = 0, {
      incProgress(0.1, detail = "leyendo secuencias")
      seqs <- read_sequences(path)
      out  <- select_orf(seqs, as.integer(input$code), input$min_len, input$keep_amb,
                         progress = function(f) incProgress(0.25, detail = paste0("marco +", f)))

      # 6b. Comparación con la versión 12 de Generate_coi12.sh (si existe)
      out$comparison <- NULL
      if (input$source == "data" && grepl("_123\\.fasta$", input$data_file)) {
        f12 <- file.path(data_dir, sub("_123\\.fasta$", "_12.fasta", input$data_file))
        if (file.exists(f12)) {
          incProgress(0.1, detail = "comparando con la versión 12")
          old <- read_sequences(f12)
          old <- setNames(as.character(old), accession(names(old)))
          new <- setNames(out$coi12, accession(names(out$coi12)))
          common <- intersect(names(new), names(old))
          out$comparison <- c(comunes = length(common),
                              iguales = sum(old[common] == new[common]))
        }
      }
      out$seqs  <- seqs
      out$table <- as.integer(input$code)
      out
    })
  })

  # 6c. Resumen
  output$summary <- renderUI({
    r <- res()$report; n <- nrow(r)
    fr <- table(factor(r$frame, 1:3))
    txt <- sprintf(
      "<ul><li><b>Secuencias leídas:</b> %s</li>
       <li><b>Conservadas:</b> %s · <b>excluidas:</b> %s</li>
       <li><b>Marco mayoritario:</b> +%d</li>
       <li><b>Fuera del marco +2:</b> %s → las que <code>Generate_coi12.sh</code> habría recortado mal</li>
       <li><b>Posibles NUMT / errores (stops en los 3 marcos):</b> %s</li></ul>",
      format(n, big.mark = ".", decimal.mark = ","), format(sum(r$kept), big.mark = ".", decimal.mark = ","),
      format(n - sum(r$kept), big.mark = ".", decimal.mark = ","), res()$majority,
      format(sum(fr[c("1", "3")]), big.mark = ".", decimal.mark = ","),
      format(sum(r$reason == "ambiguo_con_stops"), big.mark = ".", decimal.mark = ","))
    cmp <- res()$comparison
    if (!is.null(cmp))
      txt <- paste0(txt, sprintf(
        "<div class='alert alert-info'><b>Comparación con la versión 12 actual:</b><br>
         %s secuencias en común · <b>%s idénticas (%.2f %%)</b> · %s distintas</div>",
        format(cmp["comunes"], big.mark = ".", decimal.mark = ","), format(cmp["iguales"], big.mark = ".", decimal.mark = ","),
        100 * cmp["iguales"] / max(cmp["comunes"], 1),
        format(cmp["comunes"] - cmp["iguales"], big.mark = ".", decimal.mark = ",")))
    HTML(txt)
  })

  # 6d. Gráficos (R base, sin dependencias extra)
  output$p_frames <- renderPlot({
    fr <- table(factor(res()$report$frame, 1:3))
    bp <- barplot(fr, names.arg = paste0("+", 1:3), col = col_frame, border = NA,
                  ylim = c(0, max(fr, 1) * 1.25), ylab = "Nº de secuencias", las = 1)
    text(bp, fr, sprintf("%d\n(%.1f%%)", fr, 100 * fr / max(sum(fr), 1)), pos = 3, cex = 0.9)
  })
  output$p_reasons <- renderPlot({
    rs <- sort(table(res()$report$reason))
    par(mar = c(4, 13, 1, 4))
    bp <- barplot(rs, horiz = TRUE, las = 1, col = col_reason[names(rs)], border = NA,
                  xlim = c(0, max(rs) * 1.25), xlab = "Nº de secuencias", cex.names = 0.85)
    text(rs, bp, rs, pos = 4, cex = 0.9)
  })
  output$p_stops <- renderPlot({
    r <- res()$report
    boxplot(r[, c("stops_f1", "stops_f2", "stops_f3")], names = paste0("+", 1:3),
            col = col_frame, ylab = "Nº de stops internos", las = 1, outcex = 0.5)
  })

  # 6e. Desplegable del esquema y esquema
  observeEvent(res(), {
    r  <- res()$report[!is.na(res()$report$frame), ]
    ids <- r$id[order(r$frame == res()$majority)]
    updateSelectInput(session, "seq_pick", choices = head(ids, 2000), selected = ids[1])
  })
  output$scheme <- renderUI({
    req(input$seq_pick)
    r <- res()$report[res()$report$id == input$seq_pick, ][1, ]
    codon_html(r$id, as.character(res()$seqs[[r$id]]), r$frame, res()$table)
  })

  # 6f. Tabla
  output$table <- renderDT(datatable(head(res()$report, 5000), filter = "top",
                                     options = list(pageLength = 15, scrollX = TRUE)))

  # 6g. Descargas
  output$dl_fasta <- downloadHandler(
    filename = function() "coi12_orf.fasta",
    content  = function(file) writeLines(paste0(">", names(res()$coi12), "\n", res()$coi12), file))
  output$dl_report <- downloadHandler(
    filename = function() "orf_report.tsv",
    content  = function(file) write.table(res()$report, file, sep = "\t",
                                          quote = FALSE, row.names = FALSE))
}


# -----------------------------------------------------------------------------
# PASO 7. Lanzar la app
# -----------------------------------------------------------------------------
shinyApp(ui, server)