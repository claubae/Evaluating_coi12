#!/usr/bin/env python3
# =============================================================================
#  app_orf_gradio.py  —  Prototipo de interfaz web para selección de ORF y coi12
# =============================================================================
#
#  QUÉ ES
#  Una pequeña aplicación web (hecha con Gradio, igual que la del clasificador
#  de melanoma) que permite:
#     1. Subir un FASTA/CSV o pegar secuencias de COI
#     2. Elegir el código genético y las opciones del análisis
#     3. Ver qué marco de lectura se elige para cada secuencia y por qué
#     4. Ver gráficos (marcos, motivos, stops) y un esquema codón a codón
#     5. Descargar el FASTA coi12 y el informe
#
#  Reutiliza las funciones de orf.py (misma carpeta), así que los resultados
#  son exactamente los mismos que con el script de terminal.
#
#  CÓMO SE EJECUTA (desde la raíz del proyecto)
#     pip install gradio biopython pandas matplotlib
#     python3 02_Scripts/app_orf_gradio.py
#  Se abre en el navegador: http://127.0.0.1:7860
#  Para un enlace público temporal (compartirlo con el grupo): share=True
#  en demo.launch() al final del script.
#
#  DOCUMENTACIÓN
#  - Gradio, guía rápida:       https://www.gradio.app/guides/quickstart
#  - Gradio Blocks (layouts):   https://www.gradio.app/guides/blocks-and-event-listeners
#  - Abid et al. (2019) Gradio: Hassle-free sharing and testing of ML models
#    in the wild. arXiv:1906.02569. https://arxiv.org/abs/1906.02569
#  - Fundamento biológico y referencias del método: ver cabecera de orf.py
#
#  NOTA: es un PROTOTIPO. Con la base MIDORI completa (1,5 M secuencias)
#  funciona pero tarda unos minutos; para eso es mejor orf.py en terminal.
# =============================================================================


# -----------------------------------------------------------------------------
# PASO 0. Librerías y funciones de orf.py
# -----------------------------------------------------------------------------
import os
import sys
import tempfile
from collections import Counter

import gradio as gr
import matplotlib
matplotlib.use("Agg")                      # dibujar sin ventana (servidor)
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import Seq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from orf import read_sequences, count_internal_stops, longest_orf, make_coi12   # noqa: E402


# -----------------------------------------------------------------------------
# PASO 1. Constantes de la interfaz
# -----------------------------------------------------------------------------
GENETIC_CODES = {                          # https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    "5 · Mitocondrial de invertebrados (Calanoida, artrópodos)": 5,
    "2 · Mitocondrial de vertebrados": 2,
    "4 · Mitocondrial de mohos, protozoos y celentéreos": 4,
    "1 · Estándar (nuclear)": 1,
}
FRAME_COLORS = {1: "#E69F00", 2: "#0072B2", 3: "#009E73"}
REASON_COLORS = {"unico_sin_stops": "#0072B2", "empate_marco_mayoritario": "#56B4E9",
                 "empate_orf_mas_largo": "#E69F00", "ambiguo_con_stops": "#D55E00",
                 "demasiado_corta": "grey"}

EXAMPLE = (  # 2 queries reales de Calanoida + 1 desplazada 1 pb para ver el efecto
    ">JX868440_Eodiaptomus_wolterecki\nTCTTTCAAGGAATATTGCTCACGCAGGCAGTTCTGTAGATTTTGCTATTTTTTCTCTACACTTAGCCGGGGTGAGATCAATTTTAGGTGCTGTGAATTTTATTAGAACATTAGGGAATCTTCGAGCGTTCGGGATAATTTTAGATCGAATGCCTCTATTTGCTTGGGCTGTCCTAGTTACAGC\n>KC594122_Bestiolina_similis\nTCTATCTAGTAATATTGCTCATGCAGGTGCCTCAGTTGACTTTGCGATTTTTTCTTTACATTTAGCTGGTGTGAGATCAATTTTAGGAGCAGTAAATTTCATTAGTACATTAGGTAATCTTCGAGTGTTTGGAATGCTACTTGATCGTATACCATTATTTGCCTGAGCTGTTTTAATTACAGCAATTTTATTACTGCTTTCATTACCAGTGTTAGCAGGTGCAATTACTATGTTGTTAACAGATCGAAACTTAAATACATCATTTTATGATGTAGGAGGAGGAGGAGATCCTATTTTATACCAACATTTGTTTTGATTTTT\n>MW124831_Euchirella_amoena_desplazada_1pb\nAATTATCTTCTAATATTGCGCATGCAGGAGGATCTGTAGATTTTGCTATTTTCTCATTACATTTAGCAGGAGTAAGATCAATTTTAGGAGCTGTAAATTTTATTAGAACTTTAGGTAATTTACGTGTATTTGGAATGTTATTAGATCGGATACCTTTATTTGCTTGATCTGTATTAATTACAGCTATTCTTTTACTTTTATCTTTACCTGTTTTAGCTGGAGCTATTACTATACTACTAACAGATCGAAACTTAAATACAACTTTTTATGACGTTGGAGGTGGAGGGGATCCTATTTTATATCAACATTTATTT\n")


# -----------------------------------------------------------------------------
# PASO 2. Análisis (misma lógica que orf.py, pero devolviendo un DataFrame)
# -----------------------------------------------------------------------------
def analyse(records, table, min_length, keep_ambiguous):
    """records = lista de (id, secuencia). Devuelve (DataFrame informe, dict coi12)."""
    stats = []
    for _, seq in records:
        if len(seq) < min_length:
            stats.append(None)
        else:
            res = [count_internal_stops(seq, f, table) for f in range(3)]
            stats.append(([r[0] for r in res], [r[1] for r in res]))

    # Marco mayoritario entre secuencias con un único marco limpio
    uniq = Counter(st.index(0) + 1 for st, _ in filter(None, stats) if st.count(0) == 1)
    majority = uniq.most_common(1)[0][0] if uniq else 2

    rows, coi12 = [], {}
    for (name, seq), s in zip(records, stats):
        if s is None:
            rows.append(dict(id=name, length=len(seq), frame=None, reason="demasiado_corta",
                             stops_f1=None, stops_f2=None, stops_f3=None, kept=False))
            continue
        stops, prots = s
        open_f = [f + 1 for f in range(3) if stops[f] == 0]
        if len(open_f) == 1:
            frame, reason = open_f[0], "unico_sin_stops"
        elif majority in open_f:
            frame, reason = majority, "empate_marco_mayoritario"
        elif open_f:
            frame, reason = max(open_f, key=lambda f: longest_orf(prots[f - 1])), "empate_orf_mas_largo"
        else:
            frame, reason = stops.index(min(stops)) + 1, "ambiguo_con_stops"
        kept = reason != "ambiguo_con_stops" or keep_ambiguous
        rows.append(dict(id=name, length=len(seq), frame=frame, reason=reason,
                         stops_f1=stops[0], stops_f2=stops[1], stops_f3=stops[2], kept=kept))
        if kept:
            coi12[name] = make_coi12(seq, frame)
    return pd.DataFrame(rows), coi12, majority


# -----------------------------------------------------------------------------
# PASO 3. Gráficos (matplotlib -> se muestran en gr.Plot)
# -----------------------------------------------------------------------------
def plot_frames(df):
    counts = df["frame"].dropna().astype(int).value_counts().reindex([1, 2, 3], fill_value=0)
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    bars = ax.bar([f"+{f}" for f in counts.index], counts.values,
                  color=[FRAME_COLORS[f] for f in counts.index])
    total = max(counts.sum(), 1)
    ax.bar_label(bars, labels=[f"{v}\n({100 * v / total:.1f}%)" for v in counts.values])
    ax.set_ylim(0, max(counts.max(), 1) * 1.3)
    ax.set_title("Marco elegido", fontweight="bold"); ax.set_ylabel("Nº de secuencias")
    fig.tight_layout()
    return fig


def plot_reasons(df):
    c = df["reason"].value_counts()
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    bars = ax.barh(c.index[::-1], c.values[::-1], color=[REASON_COLORS[r] for r in c.index[::-1]])
    ax.bar_label(bars, padding=3)
    ax.set_xlim(0, c.max() * 1.25)
    ax.set_title("Motivo de la decisión", fontweight="bold"); ax.set_xlabel("Nº de secuencias")
    fig.tight_layout()
    return fig


def plot_stops(df):
    d = df.dropna(subset=["stops_f1"])
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    bp = ax.boxplot([d["stops_f1"], d["stops_f2"], d["stops_f3"]], patch_artist=True, widths=0.5)
    for patch, f in zip(bp["boxes"], (1, 2, 3)):
        patch.set_facecolor(FRAME_COLORS[f])
    ax.set_xticks([1, 2, 3], ["+1", "+2", "+3"])
    ax.set_title("Stops internos por marco", fontweight="bold"); ax.set_ylabel("Nº de stops")
    fig.tight_layout()
    return fig


# -----------------------------------------------------------------------------
# PASO 4. Esquema codón a codón de UNA secuencia (HTML)
# -----------------------------------------------------------------------------
def codon_html(name, seq, frame, table, n_bases=90):
    """Colorea las bases por posición en el codón y tacha la 3ª.
    Debajo muestra la traducción de los 3 marcos (los stops '*' en rojo)."""
    pos_colors = {0: "#0072B2", 1: "#56B4E9", 2: "#F0E442"}
    cells = []
    for i, b in enumerate(seq[:n_bases]):
        if i < frame - 1:
            style, deco = "background:#ddd;color:#555", "line-through"
        else:
            k = (i - (frame - 1)) % 3
            style = f"background:{pos_colors[k]};color:{'#000' if k == 2 else '#fff'}"
            deco = "line-through" if k == 2 else "none"
        gap = "margin-left:4px;" if i >= frame - 1 and (i - (frame - 1)) % 3 == 0 else ""
        cells.append(f"<wbr><span style='{style};{gap}text-decoration:{deco};text-decoration-color:red;"
                     f"text-decoration-thickness:2px;padding:2px 3px;font-family:monospace;"
                     f"font-weight:bold;border-radius:2px'>{b}</span>")
    prot_lines = []
    for f in range(3):
        _, p = count_internal_stops(seq, f, table)
        p_html = "".join(f"<span style='color:red;font-weight:bold'>*</span>" if a == "*" else a
                         for a in p[: n_bases // 3])
        mark = " ← elegido" if f + 1 == frame else ""
        prot_lines.append(f"<div style='font-family:monospace'><b>+{f + 1}</b>&nbsp; {p_html}"
                          f"<b style='color:#0072B2'>{mark}</b></div>")
    legend = ("<span style='background:#0072B2;color:#fff;padding:1px 5px'>1ª</span> "
              "<span style='background:#56B4E9;color:#fff;padding:1px 5px'>2ª</span> "
              "<span style='background:#F0E442;padding:1px 5px'>3ª (se elimina)</span> "
              "<span style='background:#ddd;padding:1px 5px'>fuera de marco</span>")
    return (f"<div style='line-height:2.1;word-break:break-all'><b>{name}</b> — marco elegido <b>+{frame}</b> "
            f"(primeras {min(n_bases, len(seq))} pb de {len(seq)})<br>{''.join(cells)}"
            f"<br><small>{legend}</small><br><br><b>Traducción de los 3 marcos</b>"
            f"{''.join(prot_lines)}<br><b>coi12:</b> "
            f"<code>{make_coi12(seq, frame)[:60]}…</code></div>")


# -----------------------------------------------------------------------------
# PASO 5. Función que conecta la interfaz con el análisis
# -----------------------------------------------------------------------------
def run(file, pasted, code_label, min_length, keep_ambiguous):
    # 5a. Recoger las secuencias (archivo subido o texto pegado)
    if file is not None:
        records = read_sequences(file if isinstance(file, str) else file.name)
    elif pasted and pasted.strip():
        tmp = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
        tmp.write(pasted if pasted.lstrip().startswith(">") else ">secuencia_1\n" + pasted)
        tmp.close()
        records = read_sequences(tmp.name)
    else:
        raise gr.Error("Sube un archivo o pega al menos una secuencia.")
    if not records:
        raise gr.Error("No se ha leído ninguna secuencia. Revisa el formato.")

    # 5b. Analizar
    table = GENETIC_CODES[code_label]
    df, coi12, majority = analyse(records, table, int(min_length), keep_ambiguous)

    # 5c. Resumen en texto
    n = len(df)
    fr = df["frame"].dropna().astype(int).value_counts()
    rs = df["reason"].value_counts()
    summary = (f"### Resultado\n"
               f"- **Secuencias leídas:** {n:,}\n"
               f"- **Conservadas:** {int(df['kept'].sum()):,} · "
               f"**excluidas:** {n - int(df['kept'].sum()):,}\n"
               f"- **Marco mayoritario:** +{majority}\n"
               f"- **Fuera del marco +2:** {int(fr.drop(2, errors='ignore').sum()):,} "
               f"→ las que `Generate_coi12.sh` (marco +2 fijo) habría recortado mal\n"
               f"- **Posibles NUMT / errores (stops en los 3 marcos):** "
               f"{int(rs.get('ambiguo_con_stops', 0)):,}")
    probs = {f"+{f}": float(fr.get(f, 0) / max(fr.sum(), 1)) for f in (1, 2, 3)}

    # 5d. Archivos para descargar
    out_dir = tempfile.mkdtemp()
    fasta_path = os.path.join(out_dir, "coi12_orf.fasta")
    with open(fasta_path, "w") as o:
        for k, s in coi12.items():
            o.write(f">{k}\n{s}\n")
    report_path = os.path.join(out_dir, "orf_report.tsv")
    df.to_csv(report_path, sep="\t", index=False)

    # 5e. Esquema de la primera secuencia "interesante" (fuera del mayoritario si la hay)
    seqs = dict(records)
    interesting = df[(df["frame"] != majority) & df["frame"].notna()]
    first = (interesting if len(interesting) else df.dropna(subset=["frame"])).iloc[0]
    html = codon_html(first["id"], seqs[first["id"]], int(first["frame"]), table)

    ids = df.dropna(subset=["frame"]).sort_values(
        "frame", key=lambda s: s != majority, ascending=False)["id"].tolist()
    state = {"seqs": seqs, "df": df, "table": table}

    return (summary, probs, plot_frames(df), plot_reasons(df), plot_stops(df),
            df.head(500), html, [fasta_path, report_path],
            gr.update(choices=ids[:2000], value=first["id"]), state)


def show_sequence(seq_id, state):
    """Redibuja el esquema cuando se elige otra secuencia en el desplegable."""
    if not state or not seq_id:
        return ""
    row = state["df"].set_index("id").loc[seq_id]
    return codon_html(seq_id, state["seqs"][seq_id], int(row["frame"]), state["table"])


# -----------------------------------------------------------------------------
# PASO 6. Construcción de la interfaz (gr.Blocks)
# -----------------------------------------------------------------------------
METHOD_MD = """
## ¿Cómo decide el marco?
1. **Traduce los 3 marcos** con el código genético elegido (por defecto, tabla 5:
   mitocondrial de invertebrados; [códigos del NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)).
2. **Cuenta los codones stop internos.** COI codifica una proteína, así que en el
   marco correcto un fragmento interno no tiene stops.
3. **Decide:**
   - un solo marco sin stops → ese;
   - varios sin stops → el **marco mayoritario** del conjunto (o el ORF más largo);
   - stops en los 3 → se excluye (posible **NUMT**, [Song et al. 2008](https://doi.org/10.1073/pnas.0803076105)).
4. **Genera coi12** quitando la 3ª posición de cada codón desde ese marco.

## Referencias
- Leray et al. (2013) *Front. Zool.* 10:34 — región Leray. [doi](https://doi.org/10.1186/1742-9994-10-34)
- Leray et al. (2022) *Environmental DNA* 4:894–907 — MIDORI2. [doi](https://doi.org/10.1002/edn3.303) · [web](https://www.reference-midori.info)
- Folmer et al. (1994) *Mol. Mar. Biol. Biotechnol.* 3:294–299 — barcode COI.
- Cock et al. (2009) *Bioinformatics* 25:1422 — Biopython. [doi](https://doi.org/10.1093/bioinformatics/btp163)
- Nugent et al. (2020) *Genome* 63:291 — paquete *coil*. [doi](https://doi.org/10.1139/gen-2019-0206) · [GitHub](https://github.com/CNuge/coil)
- Abid et al. (2019) — Gradio. [arXiv](https://arxiv.org/abs/1906.02569) · [docs](https://www.gradio.app/guides/quickstart)
"""

with gr.Blocks(title="Selección de ORF · COI") as demo:

    gr.Markdown("# 🧬 Selección del marco de lectura y coi12\n"
                "Sube tus secuencias de COI (región Leray u otra) y la app decide el marco de "
                "lectura de **cada** secuencia antes de quitar la 3ª base del codón.")
    state = gr.State()

    with gr.Tab("Analizar"):
        with gr.Row():
            # ---- Columna izquierda: entradas ----
            with gr.Column(scale=1):
                file_in = gr.File(label="Archivo FASTA o CSV (taxid,seq,accession,especie)",
                                  file_types=[".fasta", ".fa", ".fas", ".txt", ".csv"])
                pasted = gr.Textbox(label="…o pega secuencias en FASTA", lines=6,
                                    placeholder=">id\nACGT...")
                code = gr.Dropdown(list(GENETIC_CODES), value=list(GENETIC_CODES)[0],
                                   label="Código genético")
                min_len = gr.Slider(3, 300, value=9, step=1,
                                    label="Longitud mínima (pb) para evaluar")
                keep_amb = gr.Checkbox(label="Conservar secuencias con stops en los 3 marcos")
                with gr.Row():
                    gr.Button("Cargar ejemplo").click(lambda: EXAMPLE, outputs=pasted)
                    btn = gr.Button("▶ Analizar", variant="primary")

            # ---- Columna derecha: resumen ----
            with gr.Column(scale=1):
                summary = gr.Markdown()
                frame_label = gr.Label(label="Proporción de secuencias en cada marco")
                files_out = gr.File(label="Descargar resultados", file_count="multiple")

        with gr.Row():
            p1 = gr.Plot(label="Marco elegido")
            p2 = gr.Plot(label="Motivo")
            p3 = gr.Plot(label="Stops por marco")

        gr.Markdown("### Esquema codón a codón")
        seq_pick = gr.Dropdown(label="Elige una secuencia (primero las que NO están en el marco mayoritario)")
        scheme = gr.HTML()
        gr.Markdown("### Informe (primeras 500 filas)")
        table_out = gr.Dataframe(wrap=True)

    with gr.Tab("Cómo funciona"):
        gr.Markdown(METHOD_MD)

    # ---- Eventos ----
    btn.click(run, inputs=[file_in, pasted, code, min_len, keep_amb],
              outputs=[summary, frame_label, p1, p2, p3, table_out, scheme, files_out,
                       seq_pick, state])
    seq_pick.change(show_sequence, inputs=[seq_pick, state], outputs=scheme)


# -----------------------------------------------------------------------------
# PASO 7. Lanzar
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    demo.launch(theme=gr.themes.Soft())   # añade share=True para un enlace público temporal