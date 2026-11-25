#!/usr/bin/env python3
# =============================================================================
#  orf.py  —  Selección del marco de lectura (ORF) y generación de coi12
# =============================================================================
#
#  QUÉ HACE, EN UNA FRASE
#  Para cada secuencia de COI decide en qué marco de lectura está (+1, +2 o +3)
#  mirando dónde NO aparecen codones stop, y luego elimina la 3ª base de cada
#  codón usando ESE marco (en vez de asumir +2 para todas, como hace
#  Generate_coi12.sh).
#
#  USO (desde la carpeta raíz del proyecto)
#     python3 02_Scripts/orf.py data/midori_123.fasta
#     python3 02_Scripts/orf.py data/all_together_calanoida_123.fasta
#     python3 02_Scripts/orf.py data/midori_123.fasta --compare data/midori_12.fasta
#
#  SALIDAS (carpeta results/orf/)
#     <nombre>_12_orf.fasta      secuencias coi12 (sin 3ª base) con su marco
#     <nombre>_orf_report.tsv    una fila por secuencia: marco, stops, motivo
#
#  REQUISITOS
#     Python >= 3.8 y Biopython:   pip install biopython
#
#  Es la versión Python de orf_seleccion_coi12.R: mismos criterios, mismos
#  resultados (comprobado con las 545 queries de Calanoida).
#
# -----------------------------------------------------------------------------
#  FUNDAMENTO BIOLÓGICO (por qué funciona)
# -----------------------------------------------------------------------------
#  1. COI es un gen que codifica proteína. Leído en el marco correcto, un
#     fragmento INTERNO del gen (como la región Leray, ~313 pb) no puede tener
#     codones stop: si los tuviera, la proteína saldría truncada.
#     En los otros dos marcos los stops aparecen por azar (~1 de cada 20
#     codones), así que normalmente solo UN marco queda "limpio".
#
#  2. COI está en el genoma mitocondrial, que usa un código genético distinto
#     del estándar. En invertebrados (tabla 5 del NCBI) TGA no es stop sino
#     triptófano (W), y AGA/AGG codifican serina. Si se usara el código estándar
#     se contarían stops falsos. Para vertebrados se usa la tabla 2.
#       -> Códigos genéticos del NCBI:
#          https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#
#  3. Si una secuencia tiene stops en los 3 marcos puede ser un pseudogén
#     nuclear de origen mitocondrial (NUMT) o tener errores de secuenciación.
#     Estas secuencias se marcan y se excluyen.
#       -> Song et al. (2008) PNAS 105:13486-13491.
#          https://doi.org/10.1073/pnas.0803076105
#
#  4. ¿Por qué quitar la 3ª posición? Es la más variable del codón (muchos
#     cambios ahí son sinónimos) y se satura rápido, lo que puede añadir ruido
#     en comparaciones entre taxones lejanos. Esa es la hipótesis del TFG.
#       -> Xia et al. (2003) Mol. Phylogenet. Evol. 26:1-7.
#          https://doi.org/10.1016/S1055-7903(02)00326-3
#
# -----------------------------------------------------------------------------
#  REFERENCIAS DEL MARCADOR Y LA BASE DE DATOS
# -----------------------------------------------------------------------------
#  - Folmer et al. (1994) Mol. Mar. Biol. Biotechnol. 3:294-299 (barcode COI).
#  - Hebert et al. (2003) Proc. R. Soc. B 270:313-321 (DNA barcoding).
#       https://doi.org/10.1098/rspb.2002.2218
#  - Leray et al. (2013) Frontiers in Zoology 10:34 (región Leray, mlCOIintF).
#       https://doi.org/10.1186/1742-9994-10-34
#  - Leray et al. (2022) Environmental DNA 4:894-907 (MIDORI2).
#       https://doi.org/10.1002/edn3.303   |   https://www.reference-midori.info
#
#  HERRAMIENTAS Y REPOS RELACIONADOS (alternativas o para ir más allá)
#  - Biopython: Cock et al. (2009) Bioinformatics 25:1422-1423.
#       https://doi.org/10.1093/bioinformatics/btp163
#       Tutorial (traducción, SeqIO): https://biopython.org/DIST/docs/tutorial/Tutorial.html
#  - coil (paquete R que corrige marco y traduce barcodes de COI):
#       Nugent et al. (2020) Genome 63:291-305. https://doi.org/10.1139/gen-2019-0206
#       https://github.com/CNuge/coil
#  - MACSE (alineamiento de secuencias codificantes respetando el marco):
#       Ranwez et al. (2011) PLoS ONE 6:e22594. https://doi.org/10.1371/journal.pone.0022594
#
#  LIBROS DE APOYO
#  - Zvelebil & Baum (2008) Understanding Bioinformatics. Garland Science.
#    (cap. sobre genes, marcos de lectura y traducción)
#  - Antao (2022) Bioinformatics with Python Cookbook, 3ª ed. Packt.
#    (recetas con Biopython: leer FASTA, traducir, filtrar)
# =============================================================================


# -----------------------------------------------------------------------------
# PASO 0. Librerías
# -----------------------------------------------------------------------------
import argparse                    # leer opciones de la terminal
import csv                         # leer el "fasta" de Calanoida (que es un CSV)
import os                          # crear carpetas y rutas
import re                          # quitar la 3ª base con una expresión regular
import sys                         # mensajes por pantalla
from collections import Counter    # contar marcos

from Bio.Seq import Seq            # traducir ADN -> proteína
from Bio.SeqIO.FastaIO import SimpleFastaParser   # lector FASTA rápido


# -----------------------------------------------------------------------------
# PASO 1. Leer las secuencias
# -----------------------------------------------------------------------------
def read_sequences(path):
    """Devuelve una lista de (id, secuencia).

    Acepta dos formatos:
      - FASTA normal (midori_123.fasta):        >AB000119_37099 / ACGT...
      - CSV del proyecto (all_together_*.fasta): taxid,secuencia,accession,especie
    """
    with open(path) as f:
        first = f.readline()

    records = []
    if first.startswith(">"):                          # FASTA
        with open(path) as f:
            for header, seq in SimpleFastaParser(f):
                records.append((header, seq))
    else:                                              # CSV
        with open(path) as f:
            for row in csv.reader(f):
                if len(row) >= 4:
                    records.append((f"{row[2]};{row[3]}", row[1]))

    # Limpieza: mayúsculas, sin gaps y cualquier carácter raro -> N
    clean = []
    for name, seq in records:
        seq = seq.upper().replace("-", "")
        seq = re.sub(r"[^ACGTRYSWKMBDHVN]", "N", seq)
        clean.append((name, seq))
    return clean


# -----------------------------------------------------------------------------
# PASO 2. Traducir un marco y contar stops internos
# -----------------------------------------------------------------------------
def count_internal_stops(seq, frame, table):
    """frame = 0, 1 o 2 (marcos +1, +2, +3).

    1) Se saltan las primeras 'frame' bases.
    2) Se recorta el final para tener codones completos.
    3) Se traduce con el código genético indicado (los stops salen como '*').
    4) Se cuentan los '*' sin contar el último aminoácido: un stop justo al
       final puede deberse a dónde se cortó el fragmento, no a la biología.
    """
    sub = seq[frame:]
    sub = sub[: len(sub) - len(sub) % 3]
    protein = str(Seq(sub).translate(table=table))
    return protein[:-1].count("*"), protein


def longest_orf(protein):
    """Longitud (en aminoácidos) del tramo más largo sin stops."""
    return max(len(piece) for piece in protein.split("*"))


# -----------------------------------------------------------------------------
# PASO 3. Quitar la 3ª base de cada codón
# -----------------------------------------------------------------------------
def make_coi12(seq, frame):
    """frame = 1, 2 o 3. Desde el primer codón completo, el patrón '(..).'
    se queda con 2 bases y descarta la 3ª, una y otra vez.
    Ejemplo (marco +2):  A|TTA|TCT|TCT  ->  TT TC TC  ->  TTTCTC"""
    return re.sub(r"(..).", r"\1", seq[frame - 1:])


# -----------------------------------------------------------------------------
# PASO 4. Programa principal
# -----------------------------------------------------------------------------
def main():
    # 4a. Opciones
    ap = argparse.ArgumentParser(description="Selección de ORF por secuencia y coi12")
    ap.add_argument("input", help="FASTA o CSV con las secuencias 123")
    ap.add_argument("--out-dir", default="results/orf", help="carpeta de salida")
    ap.add_argument("--table", type=int, default=5,
                    help="código genético NCBI (5 = mitocondrial invertebrados)")
    ap.add_argument("--min-length", type=int, default=9,
                    help="secuencias más cortas no se evalúan")
    ap.add_argument("--keep-ambiguous", action="store_true",
                    help="conservar secuencias con stops en los 3 marcos")
    ap.add_argument("--compare", help="(opcional) FASTA/CSV 12 hecho con Generate_coi12.sh")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    base = re.sub(r"_?123\.fasta$", "", os.path.basename(a.input))

    # 4b. Leer
    seqs = read_sequences(a.input)
    print(f"Secuencias leídas: {len(seqs):,}", file=sys.stderr)

    # 4c. PRIMERA PASADA: stops en cada marco
    stats = []                                  # (stops[3], proteínas[3]) o None
    for i, (name, seq) in enumerate(seqs, 1):
        if len(seq) < a.min_length:
            stats.append(None)
        else:
            res = [count_internal_stops(seq, f, a.table) for f in range(3)]
            stats.append(([r[0] for r in res], [r[1] for r in res]))
        if i % 200_000 == 0:
            print(f"  ... {i:,} secuencias traducidas", file=sys.stderr)

    # 4d. Marco MAYORITARIO: solo cuentan las secuencias con un único marco
    #     limpio (evidencia sólida). Se usa para desempatar.
    unique_frames = Counter(st.index(0) + 1 for st, _ in filter(None, stats) if st.count(0) == 1)
    majority = unique_frames.most_common(1)[0][0]
    print(f"Marco mayoritario: +{majority}", file=sys.stderr)

    # 4e. SEGUNDA PASADA: decidir el marco de cada secuencia
    #     - 1 marco sin stops          -> ese         (unico_sin_stops)
    #     - varios sin stops           -> mayoritario (empate_marco_mayoritario)
    #                                     o ORF más largo si el mayoritario no está
    #     - stops en los 3             -> excluida    (ambiguo_con_stops, ¿NUMT?)
    rows, coi12 = [], []
    for (name, seq), s in zip(seqs, stats):
        if s is None:
            rows.append([name, len(seq), "", "demasiado_corta", "", "", "", False])
            continue
        stops, prots = s
        open_frames = [f + 1 for f in range(3) if stops[f] == 0]
        if len(open_frames) == 1:
            frame, reason = open_frames[0], "unico_sin_stops"
        elif majority in open_frames:
            frame, reason = majority, "empate_marco_mayoritario"
        elif open_frames:
            frame = max(open_frames, key=lambda f: longest_orf(prots[f - 1]))
            reason = "empate_orf_mas_largo"
        else:
            frame, reason = stops.index(min(stops)) + 1, "ambiguo_con_stops"

        kept = reason != "ambiguo_con_stops" or a.keep_ambiguous
        rows.append([name, len(seq), frame, reason, *stops, kept])
        if kept:
            coi12.append((name, make_coi12(seq, frame)))

    # 4f. Guardar FASTA coi12 (una línea por secuencia, como Generate_coi12.sh)
    fasta_out = os.path.join(a.out_dir, f"{base}_12_orf.fasta")
    with open(fasta_out, "w") as o:
        for name, s in coi12:
            o.write(f">{name}\n{s}\n")

    # 4g. Guardar informe
    report_out = os.path.join(a.out_dir, f"{base}_orf_report.tsv")
    with open(report_out, "w", newline="") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["id", "length", "frame", "reason", "stops_f1", "stops_f2", "stops_f3", "kept"])
        for r in rows:
            r[-1] = "TRUE" if r[-1] else "FALSE"
            w.writerow(r)

    # 4h. Resumen por pantalla
    frames = Counter(r[2] for r in rows if r[2] != "")
    reasons = Counter(r[3] for r in rows)
    n = len(rows)
    print("\nMarco elegido:", file=sys.stderr)
    for f in (1, 2, 3):
        print(f"  +{f}: {frames[f]:>10,}  ({100 * frames[f] / n:5.2f} %)", file=sys.stderr)
    print("Motivo:", file=sys.stderr)
    for k, v in reasons.most_common():
        print(f"  {k:<26} {v:>10,}  ({100 * v / n:5.2f} %)", file=sys.stderr)
    print(f"\nConservadas: {len(coi12):,} de {n:,}\n-> {fasta_out}\n-> {report_out}", file=sys.stderr)

    # 4i. (Opcional) Comparar con el coi12 hecho con Generate_coi12.sh (marco +2 fijo)
    if a.compare:
        # Se compara por ACCESSION (lo que va antes de '_' o ';'), porque
        # midori_123 usa ">AB000119_37099" y midori_12 solo ">AB000119".
        acc = lambda h: re.split(r"[_;\s]", h)[0]
        old = {acc(k): s for k, s in read_sequences(a.compare)}
        common = [(acc(k), s) for k, s in coi12 if acc(k) in old]
        same = sum(old[k] == s for k, s in common)
        print(f"\nComparación con {a.compare}:\n"
              f"  secuencias en común: {len(common):,}\n"
              f"  idénticas:           {same:,} ({100 * same / max(len(common), 1):.2f} %)\n"
              f"  distintas:           {len(common) - same:,}  (las que NO estaban en +2)",
              file=sys.stderr)


if __name__ == "__main__":
    main()