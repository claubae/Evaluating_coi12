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
#  SALIDAS (carpeta 08_orf_results/<nombre>/)
#     <nombre>_12_orf.fasta      secuencias coi12 (sin 3ª base) con su marco
#     <nombre>_orf_report.tsv    una fila por secuencia: marco, stops, motivo
#     estadisticas.csv           proporciones con IC 95 % (Wilson)
#     figuras_python/            7 figuras (marcos, motivos, stops, IC 95 %,
#                                longitud, identidad coi12 antiguo/nuevo, infografía)
#  Si existe la pareja *_12.fasta en la misma carpeta, compara automáticamente.
#
#  REQUISITOS
#     Python >= 3.8:   pip install biopython pandas matplotlib numpy
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
    ap.add_argument("--out-dir", default=None,
                    help="carpeta de salida (por defecto 08_orf_results/<nombre>)")
    ap.add_argument("--table", type=int, default=5,
                    help="código genético NCBI (5 = mitocondrial invertebrados)")
    ap.add_argument("--min-length", type=int, default=9,
                    help="secuencias más cortas no se evalúan")
    ap.add_argument("--keep-ambiguous", action="store_true",
                    help="conservar secuencias con stops en los 3 marcos")
    ap.add_argument("--compare", help="(opcional) FASTA/CSV 12 hecho con Generate_coi12.sh")
    a = ap.parse_args()
    base = re.sub(r"_?123\.fasta$", "", os.path.basename(a.input))
    if a.out_dir is None:
        a.out_dir = os.path.join("08_orf_results", base)
    os.makedirs(a.out_dir, exist_ok=True)
    # Si no se indica --compare, se busca la pareja *_12.fasta en la misma carpeta
    if a.compare is None and a.input.endswith("_123.fasta"):
        cand = a.input[:-len("_123.fasta")] + "_12.fasta"
        if os.path.exists(cand):
            a.compare = cand

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

    # 4i. Comparar con el coi12 hecho con Generate_coi12.sh (marco +2 fijo)
    #     Se compara por ACCESSION (lo que va antes de '_' o ';'), porque
    #     midori_123 usa ">AB000119_37099" y midori_12 solo ">AB000119".
    #     Si un accession aparece varias veces se usa la PRIMERA copia (igual que en R).
    comparison = None
    if a.compare:
        acc = lambda h: re.split(r"[_;\s]", h)[0]
        old, new = {}, {}
        for k, s in read_sequences(a.compare):
            old.setdefault(acc(k), s)
        for k, s in coi12:
            new.setdefault(acc(k), s)
        common = [k for k in new if k in old]
        diff = [k for k in common if old[k] != new[k]]
        same = len(common) - len(diff)
        print(f"\nComparación con {a.compare}:\n"
              f"  secuencias en común: {len(common):,}\n"
              f"  idénticas:           {same:,} ({100 * same / max(len(common), 1):.2f} %)\n"
              f"  distintas:           {len(diff):,}  (las que NO estaban en +2)",
              file=sys.stderr)
        comparison = dict(common=len(common), same=same, old=old, new=new, diff=diff)

    # 4j. Estadística y figuras
    make_figures(rows, coi12, comparison, base, a.out_dir)


# -----------------------------------------------------------------------------
# PASO 5. ESTADÍSTICA Y FIGURAS  (-> <out_dir>/figuras y estadisticas.csv)
# -----------------------------------------------------------------------------
def wilson(x, n, z=1.96):
    """Proporción (%) con IC 95 % de Wilson. Brown, Cai & DasGupta (2001)
    Statistical Science 16:101-133. https://doi.org/10.1214/ss/1009213286"""
    if n == 0:
        return 0.0, 0.0, 0.0
    p = x / n
    d = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / d
    m = z * ((p * (1 - p) / n + z ** 2 / (4 * n ** 2)) ** 0.5) / d
    return 100 * p, 100 * max(0, c - m), 100 * min(1, c + m)


def make_figures(rows, coi12, comparison, base, out_dir):
    import matplotlib
    matplotlib.use("Agg")                     # dibujar sin ventana
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    fig_dir = os.path.join(out_dir, "figuras_python")   # las de R van en figuras/
    os.makedirs(fig_dir, exist_ok=True)
    C_FR = {1: "#E69F00", 2: "#0072B2", 3: "#009E73"}
    C_RS = {"unico_sin_stops": "#0072B2", "empate_marco_mayoritario": "#56B4E9",
            "empate_orf_mas_largo": "#E69F00", "ambiguo_con_stops": "#D55E00",
            "demasiado_corta": "grey"}
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "font.size": 11})
    fmt = lambda x: f"{x:,}".replace(",", ".")

    df = pd.DataFrame(rows, columns=["id", "length", "frame", "reason",
                                     "stops_f1", "stops_f2", "stops_f3", "kept"])
    df["kept"] = df["kept"].isin([True, "TRUE"])
    ev = df[df["reason"] != "demasiado_corta"].copy()
    ev["frame"] = ev["frame"].astype(int)
    for c in ("stops_f1", "stops_f2", "stops_f3"):
        ev[c] = ev[c].astype(int)
    n_ev = len(ev)

    # 5a. Tabla de estadísticas con IC 95 %
    cats = {
        "Marco +2 (el asumido)": (ev.frame == 2) & ev.kept,
        "Marco +1 (corregidas)": (ev.frame == 1) & ev.kept,
        "Marco +3 (corregidas)": (ev.frame == 3) & ev.kept,
        "Fuera de +2 (total corregidas)": (ev.frame != 2) & ev.kept,
        "Excluidas: stops en los 3 marcos": ev.reason == "ambiguo_con_stops",
        "Resueltas por desempate": ev.reason.isin(["empate_marco_mayoritario", "empate_orf_mas_largo"]),
    }
    est = []
    for k, mask in cats.items():
        x = int(mask.sum())
        p, lo, hi = wilson(x, n_ev)
        est.append(dict(categoria=k, n=x, total=n_ev, pct=p, ic_inf=lo, ic_sup=hi))
    if comparison:
        p, lo, hi = wilson(comparison["same"], comparison["common"])
        est.append(dict(categoria="Idénticas a la versión 12 actual", n=comparison["same"],
                        total=comparison["common"], pct=p, ic_inf=lo, ic_sup=hi))
    est = pd.DataFrame(est)
    est.to_csv(os.path.join(out_dir, "estadisticas.csv"), index=False)
    print("\n" + est.to_string(index=False), file=sys.stderr)

    # 5b. FIGURA 01: marco elegido
    fr = ev.frame.value_counts().reindex([1, 2, 3], fill_value=0)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    b = ax.bar([f"+{f}" for f in fr.index], fr.values, color=[C_FR[f] for f in fr.index])
    ax.bar_label(b, labels=[f"{fmt(v)}\n({100 * v / max(n_ev, 1):.3g} %)" for v in fr.values])
    ax.set_ylim(0, fr.max() * 1.25); ax.set_ylabel("Nº de secuencias")
    ax.set_title("Marco de lectura elegido")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "01_barras_marco.png"), dpi=300); plt.close(fig)

    # 5c. FIGURA 02: motivo de la decisión
    rs = df.reason.value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(7.5, 4))
    b = ax.barh(rs.index, rs.values, color=[C_RS[r] for r in rs.index])
    ax.bar_label(b, labels=[f"{fmt(v)} ({100 * v / len(df):.3g} %)" for v in rs.values], padding=3)
    ax.set_xlim(0, rs.max() * 1.35); ax.set_xlabel("Nº de secuencias")
    ax.set_title("Motivo de la decisión")
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "02_barras_motivo.png"), dpi=300); plt.close(fig)

    # 5d. FIGURA 03: stops internos por marco (muestra de hasta 50.000 para que sea rápido)
    smp = ev.sample(min(n_ev, 50_000), random_state=1)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    bp = ax.boxplot([smp.stops_f1, smp.stops_f2, smp.stops_f3], widths=0.45,
                    patch_artist=True, flierprops=dict(markersize=2),
                    showmeans=True, meanprops=dict(marker="D", mfc="white", mec="black"))
    for patch, f in zip(bp["boxes"], (1, 2, 3)):
        patch.set_facecolor(C_FR[f])
    ax.set_xticks([1, 2, 3], ["+1", "+2", "+3"]); ax.set_xlabel("Marco")
    ax.set_ylabel("Nº de stops internos")
    ax.set_title("Codones stop internos en cada marco\n(◇ = media; el marco correcto queda en 0)", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "03_stops_por_marco.png"), dpi=300); plt.close(fig)

    # 5e. FIGURA 04: proporciones con IC 95 % (escala log)
    e = est[(est.n > 0) & ~est.categoria.str.startswith("Idénticas")].iloc[::-1]
    col = ["#0072B2" if "asumido" in c else "#D55E00" if "Excluidas" in c else "#E69F00"
           for c in e.categoria]
    fig, ax = plt.subplots(figsize=(10, 5.5))
    y = np.arange(len(e))
    for yi, (p, lo, hi, c) in enumerate(zip(e.pct, e.ic_inf.clip(lower=1e-6), e.ic_sup, col)):
        ax.errorbar(p, yi, xerr=[[p - lo], [hi - p]], fmt="none", ecolor=c, capsize=5, lw=1.5)
    ax.scatter(e.pct, y, c=col, s=70, zorder=3)
    for yi, (n, p, c) in enumerate(zip(e.n, e.pct, col)):
        ptxt = (f"{p:.2g}" if p < 99 else f"{p:.6g}").replace(".", ",")
        ax.annotate(f"{fmt(int(n))} sec. · {ptxt} %", (p, yi), xytext=(0, 9),
                    textcoords="offset points", ha="center", color=c, fontsize=9)
    ax.set_xscale("log"); ax.set_yticks(y, e.categoria)
    ax.set_xlabel("% de secuencias (escala logarítmica)")
    ax.set_title(f"¿Cuántas secuencias NO están en el marco +2?\n{fmt(n_ev)} evaluadas · punto = proporción, barra = IC 95 % (Wilson)", fontsize=12)
    ax.margins(x=0.15, y=0.12)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "04_proporciones_IC95.png"), dpi=300); plt.close(fig)

    # 5f. FIGURA 05: fuera de +2 según la longitud
    bins = [0, 200, 300, 310, 312, 313, 314, 316, 330, np.inf]
    labs = ["<200", "200-299", "300-309", "310-311", "312", "313", "314-315", "316-329", "≥330"]
    k = ev[ev.kept].copy()
    k["rango"] = pd.cut(k.length, bins, labels=labs, right=False)
    g = k.groupby("rango", observed=True).agg(n=("frame", "size"),
                                              fuera=("frame", lambda s: int((s != 2).sum())))
    g = g[g.n >= 30]
    w = [wilson(f, n) for f, n in zip(g.fuera, g.n)]
    g["pct"], g["ic_inf"], g["ic_sup"] = zip(*w) if len(w) else ([], [], [])
    g.to_csv(os.path.join(out_dir, "fuera_de_marco_por_longitud.csv"))
    fig, ax = plt.subplots(figsize=(9, 5))
    b = ax.bar(g.index.astype(str), g.pct, color="#E69F00", alpha=0.85,
               yerr=[g.pct - g.ic_inf, g.ic_sup - g.pct], capsize=4)
    for xi, (f, n, hi) in enumerate(zip(g.fuera, g.n, g.ic_sup)):
        ax.text(xi, hi, f"{fmt(int(f))}\nde {fmt(int(n))}", ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Longitud de la secuencia (pb)"); ax.set_ylabel("% fuera del marco +2")
    ax.set_title("Secuencias fuera del marco +2 según su longitud (IC 95 %)")
    ax.margins(y=0.15)
    fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "05_fuera_de_marco_por_longitud.png"), dpi=300); plt.close(fig)

    # 5g. FIGURA 06: identidad coi12 antiguo vs corregido
    if comparison and comparison["diff"]:
        def ident(x, y):
            L = min(len(x), len(y))
            return 100 * sum(c1 == c2 for c1, c2 in zip(x[:L], y[:L])) / L if L else np.nan
        old, new = comparison["old"], comparison["new"]
        frame_of = {re.split(r"[_;\s]", r[0])[0]: r[2] for r in rows}
        idt = pd.DataFrame([dict(accession=k, identidad=ident(old[k], new[k]),
                                 marco=frame_of.get(k)) for k in comparison["diff"]])
        idt.to_csv(os.path.join(out_dir, "identidad_coi12_antiguo_vs_corregido.csv"), index=False)
        fig, ax = plt.subplots(figsize=(8, 4.8))
        groups = [idt[idt.marco == f].identidad for f in (1, 3)]
        ax.hist(groups, bins=np.arange(0, 102.5, 2.5), stacked=True,
                color=[C_FR[1], C_FR[3]], label=["Marco real +1", "Marco real +3"], edgecolor="white")
        med = idt.identidad.median()
        ax.axvline(med, ls="--", c="k"); ax.text(med, ax.get_ylim()[1] * 0.95, f" mediana = {med:.1f} %")
        ax.set_xlim(0, 100); ax.legend()
        ax.set_xlabel("Identidad posición a posición (%)"); ax.set_ylabel("Nº de secuencias")
        ax.set_title(f"Parecido entre el coi12 antiguo y el corregido ({len(idt)} secuencias)")
        fig.tight_layout(); fig.savefig(os.path.join(fig_dir, "06_identidad_coi12_antiguo_vs_corregido.png"), dpi=300); plt.close(fig)

    # 5h. FIGURA 07: infografía resumen
    pct2 = est.loc[est.categoria == "Marco +2 (el asumido)", "pct"].iloc[0]
    cards = [
        (fmt(len(df)), "secuencias analizadas", "#4D4D4D"),
        (f"{pct2:.3f} %".replace(".", ","), "en el marco +2\n(la suposición del TFG)", "#0072B2"),
        (fmt(int(cats["Fuera de +2 (total corregidas)"].sum())), "corregidas\n(estaban en +1 o +3)", "#E69F00"),
        (fmt(int(cats["Excluidas: stops en los 3 marcos"].sum())), "excluidas\n(stops en los 3 marcos, ¿NUMT?)", "#D55E00"),
        (fmt(int(cats["Resueltas por desempate"].sum())), "resueltas por desempate\n(varios marcos sin stops)", "#56B4E9"),
        ((f"{100 * comparison['same'] / max(comparison['common'], 1):.2f} %".replace(".", ",")
          if comparison else "—"), "idénticas al coi12 de\nGenerate_coi12.sh", "#009E73"),
    ]
    fig, ax = plt.subplots(figsize=(12, 6)); ax.axis("off")
    for i, (num, txt, c) in enumerate(cards):
        x, y = (i % 3) * 1.0, 1 - (i // 3) * 1.0
        ax.add_patch(plt.Rectangle((x + 0.03, y + 0.05), 0.94, 0.9, fc=c, alpha=0.1, ec=c, lw=2))
        ax.text(x + 0.5, y + 0.6, num, ha="center", va="center", fontsize=30, fontweight="bold", color=c)
        ax.text(x + 0.5, y + 0.28, txt, ha="center", va="center", fontsize=12, color="#404040")
    ax.set_xlim(0, 3); ax.set_ylim(0, 2)
    fig.suptitle(f"Selección del marco de lectura · {base}", fontsize=18, fontweight="bold")
    fig.text(0.5, 0.9, "Código mitocondrial de invertebrados (tabla 5) · marco = el que no tiene codones stop internos",
             ha="center", color="grey")
    fig.savefig(os.path.join(fig_dir, "07_infografia_resumen.png"), dpi=300, bbox_inches="tight"); plt.close(fig)

    print(f"Figuras y estadísticas en: {out_dir}", file=sys.stderr)

if __name__ == "__main__":
    main()