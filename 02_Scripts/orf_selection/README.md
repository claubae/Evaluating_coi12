# Selección del marco de lectura (ORF) y coi12

Aportación al proyecto **Evaluating_coi12**. Contiene scripts en R y Python, y dos interfaces web, que **deciden el marco de lectura de cada secuencia de COI** antes de quitar la 3ª base del codón. Así no hace falta asumir que todas están en el marco +2, como hace `Generate_coi12.sh`.

> **Resultado principal:** en MIDORI (1.518.144 secuencias), el **99,99 %** está en el marco +2, así que la suposición del TFG era correcta en casi todos los casos. Hay **96 secuencias en otro marco** (64 en +1 y 32 en +3), que `Generate_coi12.sh` recortaba mal y estos scripts corrigen, y **27 posibles NUMT o errores**, que quedan excluidos. En las 545 queries de Calanoida, el 100 % está en +2.

---

## 1. ¿Por qué hace falta?

Para generar coi12 se elimina la 3ª posición de cada codón. Si una secuencia está **desplazada** respecto al marco asumido (por indels, por un recorte distinto en cutadapt o por errores de secuenciación), se eliminaría la 1ª o la 2ª base en lugar de la 3ª, y la secuencia coi12 de esa muestra dejaría de tener sentido biológico.

```
Marco +2 (correcto):   A | TTA | TCT | TCT ...   ->  TT TC TC  ->  TTTCTC
Si estuviera en +3:    AT | TAT | CTT | CT ...    ->  con +2 se quitaría la base equivocada
```

## 2. Cómo decide el marco

| Paso | Criterio | Fundamento |
|---|---|---|
| 1 | Traduce los 3 marcos con el **código mitocondrial de invertebrados (tabla 5)** | En la mitocondria TGA no es stop sino triptófano. [Códigos genéticos del NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) |
| 2 | Cuenta los **codones stop internos** (ignora el último codón) | La región Leray es un fragmento interno de un gen que codifica proteína: en el marco correcto no puede tener stops |
| 3 | **Un solo marco sin stops** → ese | Es el caso del 97,4 % de MIDORI |
| 4 | **Varios marcos sin stops** → el **marco mayoritario** del conjunto; si no está entre ellos, el ORF más largo | Todas las secuencias se recortaron con el mismo primer |
| 5 | **Stops en los 3 marcos** → se excluye | Posible pseudogén nuclear (NUMT) o error. [Song et al. 2008](https://doi.org/10.1073/pnas.0803076105) |
| 6 | Quita la 3ª base de cada codón **desde el marco elegido** | Igual que `Generate_coi12.sh`, pero con el marco de cada secuencia |

## 3. Archivos

| Archivo | Lenguaje | Para qué |
|---|---|---|
| `orf_seleccion_coi12.R` | R | Script principal. Lee de `data/`, compara con la versión 12 actual y genera figuras |
| `orf.py` | Python | Lo mismo desde la terminal (mismos resultados) |
| `app_orf_shiny.R` | R (Shiny) | Interfaz web: elegir archivo de `data/`, ver gráficos y el esquema codón a codón, descargar |
| `app_orf_gradio.py` | Python (Gradio) | La misma interfaz en Python |
| `requirements.txt` | — | Dependencias de los scripts de Python |

## 4. Datos de entrada (carpeta `data/`)

| Archivo | Formato | Contenido |
|---|---|---|
| `all_together_calanoida_123.fasta` | CSV: `taxid,secuencia,accession,especie` | Queries de Calanoida (región Leray, sin primers) |
| `all_together_calanoida_12.fasta` | CSV | Su versión coi12 actual (para comparar) |
| `midori_123.fasta` | FASTA (`>accession_taxid`) | Base MIDORI recortada con cutadapt, **antes** de `Generate_coi12.sh` |
| `midori_12.fasta` | FASTA (`>accession`) | Su versión coi12 actual (para comparar) |

Los scripts reconocen automáticamente si el archivo es FASTA o CSV.
⚠️ `data/` **no se sube a GitHub**: los archivos de MIDORI ocupan ~800 MB, por encima del límite de 100 MB por archivo.

## 5. Uso

### R (recomendado)
```r
# 1. Abrir Evaluating_coi12.Rproj en RStudio
# 2. Instalar paquetes (solo la primera vez)
install.packages(c("here", "readr", "dplyr", "ggplot2", "tidyr", "patchwork"))
install.packages("BiocManager"); BiocManager::install("Biostrings")
# 3. Abrir 02_Scripts/orf_seleccion_coi12.R y en el PASO 1 elegir:
archivo <- "all_together_calanoida_123.fasta"   # o "midori_123.fasta"
# 4. Ejecutar todo (Ctrl+Shift+S) o bloque a bloque (Ctrl+Enter)
```

### Python
```bash
pip install -r 02_Scripts/requirements.txt
python3 02_Scripts/orf.py data/all_together_calanoida_123.fasta --compare data/all_together_calanoida_12.fasta
python3 02_Scripts/orf.py data/midori_123.fasta --compare data/midori_12.fasta
```

### Interfaces web
```r
# R: abrir 02_Scripts/app_orf_shiny.R en RStudio y pulsar "Run App"
install.packages(c("shiny", "bslib", "DT"))
```
```bash
# Python: se abre en http://127.0.0.1:7860
python3 02_Scripts/app_orf_gradio.py
```

## 6. Salidas (`results/orf/`)

| Archivo | Contenido |
|---|---|
| `<nombre>_12_orf.fasta` | Secuencias coi12, cada una con su marco, una línea por secuencia |
| `<nombre>_orf_report.tsv` | Una fila por secuencia: `id, length, frame, reason, stops_f1, stops_f2, stops_f3, kept` |
| `<nombre>_distintas_vs_12.tsv` | Secuencias que cambian respecto a `Generate_coi12.sh` |
| `figuras/<nombre>/` | 8 figuras: esquema coi12, árbol de decisión, barras, stops, mapa de calor, longitudes y panel resumen |

Valores de la columna `reason`: `unico_sin_stops`, `empate_marco_mayoritario`, `empate_orf_mas_largo`, `ambiguo_con_stops` (excluida, posible NUMT) y `demasiado_corta` (menos de 9 pb).

## 7. Resultados con los datos del proyecto

| Conjunto | Secuencias | Marco +2 | Otros marcos | Excluidas | Idénticas a la versión 12 actual | Tiempo (R) |
|---|---|---|---|---|---|---|
| Calanoida (queries) | 545 | 544 (100 %) | 0 | 1 (6 pb) | 544 / 544 (100 %) | segundos |
| MIDORI (referencia) | 1.518.144 | 1.518.042 (99,99 %) | 96 | 27 (posibles NUMT) | 1.518.021 / 1.518.117 (99,99 %) | ~2,5 min |

**Interpretación:** la suposición del marco +2 está justificada empíricamente. La selección por secuencia corrige 96 referencias mal recortadas y elimina 27 secuencias dudosas, lo que da una base coi12 algo más limpia.

## 8. Referencias

**Marcador y base de datos**
- Folmer O. et al. (1994). DNA primers for amplification of mitochondrial cytochrome c oxidase subunit I from diverse metazoan invertebrates. *Mol. Mar. Biol. Biotechnol.* 3:294–299.
- Hebert P.D.N. et al. (2003). Biological identifications through DNA barcodes. *Proc. R. Soc. B* 270:313–321. https://doi.org/10.1098/rspb.2002.2218
- Leray M. et al. (2013). A new versatile primer set targeting a short fragment of the mitochondrial COI region. *Front. Zool.* 10:34. https://doi.org/10.1186/1742-9994-10-34
- Leray M. et al. (2022). MIDORI2: a collection of quality controlled, preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic mitochondrial sequences. *Environmental DNA* 4:894–907. https://doi.org/10.1002/edn3.303 · https://www.reference-midori.info

**Biología del método**
- Song H. et al. (2008). Many species in one: DNA barcoding overestimates the number of species when nuclear mitochondrial pseudogenes are coamplified. *PNAS* 105:13486–13491. https://doi.org/10.1073/pnas.0803076105
- Xia X. et al. (2003). An index of substitution saturation and its application. *Mol. Phylogenet. Evol.* 26:1–7. https://doi.org/10.1016/S1055-7903(02)00326-3

**Herramientas**
- Biostrings (Bioconductor): https://bioconductor.org/packages/Biostrings
- Cock P.J.A. et al. (2009). Biopython. *Bioinformatics* 25:1422–1423. https://doi.org/10.1093/bioinformatics/btp163
- Nugent C.M. et al. (2020). coil: an R package for cytochrome c oxidase I (COI) DNA barcode data cleaning, translation, and error evaluation. *Genome* 63:291–305. https://doi.org/10.1139/gen-2019-0206 · https://github.com/CNuge/coil
- Ranwez V. et al. (2011). MACSE: Multiple Alignment of Coding SEquences. *PLoS ONE* 6:e22594. https://doi.org/10.1371/journal.pone.0022594
- Shiny: https://shiny.posit.co · Wickham H. (2021) *Mastering Shiny*, O'Reilly (libro gratis: https://mastering-shiny.org)
- Gradio: Abid A. et al. (2019) arXiv:1906.02569 · https://www.gradio.app

**Libros de apoyo**
- Zvelebil M. & Baum J.O. (2008). *Understanding Bioinformatics*. Garland Science.
- Antao T. (2022). *Bioinformatics with Python Cookbook*, 3ª ed. Packt.