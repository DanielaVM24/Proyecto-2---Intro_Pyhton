# -*- coding: utf-8 -*-
"""
GC% en 16S rRNA — Dashboard interactivo (Proyecto 2)

Acciones que se esperan del agente

El agente debe ser capaz de realizar la siguiente tarea:

Escribir un párrafo en el que se describen los resultados obtenidos en el análisis y se responde a la pregunta:
¿Cuál organismo tiene mayor cantidad de GC y cuál menor?

Agregar las líneas de código necesarias para realizar esta funcionalidad e integrar el resultado en el tablero.

Seleccione cinco microorganismos adicionales para 16S rRNA y agréguelos al conjunto de prueba del programa, con el fin de someterlo a mayor exigencia.

Plantear posibles fragilidades del código y proponer estrategias para corregir o mitigar dichas vulnerabilidades.

Desarrollar una presentación en la que se explican los retos técnicos, la lógica de solución de los problemas y los puntos de mejora del programa.

Documentar el código en la medida de lo posible en GitHub.

-------------------------------------------------------------
Pestañas del dashboard:
- Carga: subir FASTA, elegir umbrales de calidad, métricas y cobertura.
- Tabla: tabla por entrada (EDITABLE con casillas) y resumen por organismo; descarga CSV.
- Gráfico: barras comparativas; descarga PNG.
- Conclusión: párrafo automático (mayor/menor %GC) + descarga TXT.
- Guion: puntos para presentar (retos, lógica y mejoras).
- Fragilidades: riesgos y mitigaciones.

Ejecución:
  py -3.12 -m pip install streamlit pandas matplotlib
  py -3.12 -m streamlit run streamlit_app.py
"""

import io
import re
from typing import List, Dict, Tuple

import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

# --------- Ajustes de página ---------
st.set_page_config(page_title="GC% 16S rRNA - Dashboard", layout="wide")
plt.rcParams["font.family"] = "Times New Roman"

st.title("GC% en 16S rRNA — Dashboard interactivo")
st.caption("Sube archivos FASTA, ajusta los umbrales de calidad, compara %GC y genera una conclusión automática.")

# --------- Conjuntos de organismos (para cobertura) ---------
BASE_9 = {
    "Mycobacterium tuberculosis", "Corynebacterium diphtheriae", "Bacillus subtilis",
    "Escherichia coli", "Salmonella enterica", "Pseudomonas aeruginosa",
    "Clostridium botulinum", "Borrelia burgdorferi", "Mycoplasma genitalium"
}
EXTRA_5 = {
    "Staphylococcus aureus", "Vibrio cholerae", "Helicobacter pylori",
    "Streptococcus pneumoniae", "Listeria monocytogenes"
}

# --------- Accesiones (opcional) ---------
ACCESIONES: Dict[str, str] = {
    "Mycobacterium tuberculosis": "NR_102810.2",
    "Corynebacterium diphtheriae": "NR_037079.1",
    "Bacillus subtilis": "NR_112116.2",
    "Escherichia coli": "NR_024570.1",
    "Salmonella enterica": "NR_074910.1",
    "Pseudomonas aeruginosa": "NR_074828.1",
    "Clostridium botulinum": "NR_029157.1",
    "Borrelia burgdorferi": "NR_044732.2",
    "Mycoplasma genitalium": "NR_026155.1",
    # extra sugeridos
    "Staphylococcus aureus": "NR_113956.1",
    "Vibrio cholerae": "NR_118568.1",
    "Helicobacter pylori": "NR_028911.1",
    "Streptococcus pneumoniae": "NR_114926.1",
    "Listeria monocytogenes": "NR_074996.1",
}

# --------- Utilidades ---------
def read_fasta_text(text: str) -> str:
    """Lee FASTA, elimina encabezados y no-ACGTN; devuelve secuencia concatenada en mayúsculas."""
    seq_lines: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            continue
        seq_lines.append(re.sub(r"[^ACGTNacgtn]", "", line.strip()))
        # Solo se mantienen A/C/G/T/N; se ignoran espacios y otros símbolos
    return "".join(seq_lines).upper()

def gc_n_content(seq: str) -> Tuple[float, float, int, int, int, int, int, int]:
    """
    Calcula:
    - %GC sobre A/C/G/T
    - %N sobre (A+C+G+T+N)
    Devuelve: (%GC, %N, A, C, G, T, N_total, length_acgt)
    """
    a = seq.count("A"); c = seq.count("C"); g = seq.count("G"); t = seq.count("T"); n = seq.count("N")
    acgt = a + c + g + t
    total = acgt + n
    gc_pct = float("nan") if acgt == 0 else (g + c) * 100.0 / acgt
    n_pct  = float("nan") if total == 0 else n * 100.0 / total
    return gc_pct, n_pct, a, c, g, t, n, acgt

def infer_organism_from_filename(name: str) -> str:
    """Usa el nombre de archivo como organismo (admite _ y -)."""
    base = name.rsplit(".", 1)[0]
    base = re.sub(r"[_\-]+", " ", base).strip()
    toks = base.split()
    if len(toks) >= 2:
        return toks[0].capitalize() + " " + toks[1].lower()
    return base

def build_df_from_uploads(files, len_min: int, len_max: int, amb_threshold: float) -> pd.DataFrame:
    """Tabla por ENTRADA aplicando flags con los umbrales dados."""
    rows = []
    for uf in files:
        text = uf.getvalue().decode("utf-8", errors="ignore")
        seq = read_fasta_text(text)
        gc_pct, n_pct, a, c, g, t, n, length = gc_n_content(seq)
        organismo = infer_organism_from_filename(uf.name)
        acceso = ACCESIONES.get(organismo, "")
        flag_len = not (len_min <= length <= len_max)
        flag_amb = (n_pct >= amb_threshold) if pd.notna(n_pct) else True  # pon 0.0 para “cualquier N > 0”
        rows.append({
            "Organismo": organismo,
            "Acceso_NCBI": acceso,
            "Archivo": uf.name,
            "Longitud_bases_ACGT": length,
            "A": a, "C": c, "G": g, "T": t, "N": n,
            "%GC": round(gc_pct, 3) if pd.notna(gc_pct) else gc_pct,
            "%N": round(n_pct, 3) if pd.notna(n_pct) else n_pct,
            "Flag_Longitud": flag_len,
            "Flag_Ambiguedad": flag_amb,
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("%GC", ascending=False).reset_index(drop=True)
    return df

def summarize_by_organism(df_entries: pd.DataFrame) -> pd.DataFrame:
    """Promedia por organismo cuando hay varias entradas/archivos del mismo."""
    if df_entries.empty:
        return pd.DataFrame(columns=["Organismo","n_entradas","Longitud_total_ACGT","GC_promedio","N_promedio",
                                     "Flag_Longitud_any","Flag_Ambiguedad_any"])
    df_sum = (
        df_entries.groupby("Organismo", as_index=False)
        .agg({
            "Archivo": "count",
            "Longitud_bases_ACGT": "sum",
            "%GC": "mean",
            "%N": "mean",
            "Flag_Longitud": "any",
            "Flag_Ambiguedad": "any",
        })
        .rename(columns={
            "Archivo": "n_entradas",
            "Longitud_bases_ACGT": "Longitud_total_ACGT",
            "%GC": "GC_promedio",
            "%N": "N_promedio",
            "Flag_Longitud": "Flag_Longitud_any",
            "Flag_Ambiguedad": "Flag_Ambiguedad_any",
        })
        .sort_values("GC_promedio", ascending=False, ignore_index=True)
    )
    df_sum["GC_promedio"] = df_sum["GC_promedio"].round(3)
    df_sum["N_promedio"]  = df_sum["N_promedio"].round(3)
    return df_sum

def plot_bar(df_plot: pd.DataFrame, title: str = "Contenido de GC en 16S rRNA"):
    """Gráfico de barras con %GC promedio por organismo (o %GC simple si no hay resumen)."""
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(df_plot["Organismo"], df_plot["GC_promedio"] if "GC_promedio" in df_plot.columns else df_plot["%GC"])
    ax.set_title(title)
    ax.set_xlabel("Organismo")
    ax.set_ylabel("% GC")
    ax.tick_params(axis="x", rotation=45)
    vals = df_plot["GC_promedio"] if "GC_promedio" in df_plot.columns else df_plot["%GC"]
    for i, v in enumerate(vals):
        try:
            ax.text(i, float(v) + 0.5, f"{float(v):.1f}%", ha="center", va="bottom", fontsize=8)
        except Exception:
            pass
    plt.tight_layout()
    return fig

def auto_conclusion(df_sum: pd.DataFrame) -> str:
    """Genera el párrafo automático que responde: quién tiene mayor %GC y quién menor."""
    if df_sum.empty or df_sum["GC_promedio"].isna().all():
        return "Sube datos válidos para generar una conclusión."
    max_row = df_sum.loc[df_sum["GC_promedio"].idxmax()]
    min_row = df_sum.loc[df_sum["GC_promedio"].idxmin()]
    delta = float(max_row["GC_promedio"]) - float(min_row["GC_promedio"])
    mean_gc = float(df_sum["GC_promedio"].mean())
    return (
        f"**Mayor %GC:** *{max_row['Organismo']}* ({max_row['GC_promedio']:.2f}%).  \n"
        f"**Menor %GC:** *{min_row['Organismo']}* ({min_row['GC_promedio']:.2f}%).  \n"
        f"**Rango:** ~{delta:.2f} puntos porcentuales (promedio general: {mean_gc:.2f}%).\n\n"
        "Desde la perspectiva biológica, un **%GC más alto** suele asociarse con **mayor estabilidad** del ADN/ARN "
        "(triple enlace GC), posibles **sesgos de codón** y, en algunos grupos, **tolerancia ambiental**. "
        "Un **%GC más bajo** aparece a menudo en genomas **reducidos** y **dependencia del huésped** "
        "(p. ej., *Mycoplasma*, *Borrelia*), indicando **vías simplificadas** y nichos específicos."
    )

# --------- UI: pestañas ---------
tab_carga, tab_tabla, tab_grafico, tab_conclusion, tab_guion, tab_frag = st.tabs(
    ["📥 Carga", "📊 Tabla", "📈 Gráfico", "🧠 Conclusión", "🗂️ Guion", "🧯 Fragilidades"]
)

with tab_carga:
    st.subheader("1) Cargar archivos FASTA")
    st.caption("Formatos: .fasta, .fa, .fna (uno o varios por organismo).")
    files = st.file_uploader(
        "Arrastra o selecciona uno o varios archivos",
        type=["fasta", "fa", "fna"],
        accept_multiple_files=True
    )

    # --- Umbrales ajustables ---
    c1, c2, c3 = st.columns(3)
    with c1:
        len_min = st.number_input("Longitud mínima 16S (bp)", min_value=200, max_value=3000, value=1200, step=10)
    with c2:
        len_max = st.number_input("Longitud máxima 16S (bp)", min_value=200, max_value=3000, value=1700, step=10)
    with c3:
        amb_threshold = st.number_input(
            "%N máximo permitido", min_value=0.0, max_value=100.0, value=5.0, step=0.1,
            help="Pon 0.0 si quieres marcar cualquier N>0 como problema."
        )

    st.markdown(
        f"**Organismos base (9):** {', '.join(sorted(BASE_9))}  \n"
        f"**Adicionales sugeridos (5):** {', '.join(sorted(EXTRA_5))}"
    )

    if files:
        df_entries = build_df_from_uploads(files, len_min=len_min, len_max=len_max, amb_threshold=amb_threshold)
        df_summary = summarize_by_organism(df_entries)

        # Métricas rápidas
        c1, c2, c3, c4 = st.columns(4)
        with c1: st.metric("Entradas FASTA", len(df_entries))
        with c2: st.metric("Organismos únicos", df_summary["Organismo"].nunique())
        with c3: st.metric("Promedio %GC", f"{df_summary['GC_promedio'].mean():.2f}%")
        with c4: st.metric(
            "Rango %GC",
            f"{(df_summary['GC_promedio'].max() - df_summary['GC_promedio'].min()):.2f} pp"
        )

        # Resumen de calidad
        bad_len = df_entries[df_entries["Flag_Longitud"]]
        bad_amb = df_entries[df_entries["Flag_Ambiguedad"]]
        st.info(
            f"Entradas: {len(df_entries)}  |  Fuera de longitud: {len(bad_len)}  |  "
            f"Alta ambigüedad (≥ {amb_threshold}% N): {len(bad_amb)}"
        )
        if len(bad_len) or len(bad_amb):
            st.warning("Hay entradas con problemas (se pueden ver y editar en la pestaña **Tabla**).")
        else:
            st.success("Todas las entradas pasan los criterios de calidad actuales.")

        # Cobertura 9+5
        orgs = set(df_summary["Organismo"])
        st.info(f"Cobertura — Base: {len(BASE_9 & orgs)}/9 · Adicionales: {len(EXTRA_5 & orgs)}/5")

        # Guardar en sesión
        st.session_state["df_entries"] = df_entries
        st.session_state["df_summary"] = df_summary
        st.success("Datos cargados. Revisa las pestañas **Tabla**, **Gráfico** y **Conclusión**.")
    else:
        st.info("Aún no has subido archivos.")

with tab_tabla:
    st.subheader("2) Tablas")
    df_entries = st.session_state.get("df_entries", pd.DataFrame())
    df_summary = st.session_state.get("df_summary", pd.DataFrame())

    if df_entries.empty:
        st.warning("Sube archivos en la pestaña **Carga**.")
    else:
        st.markdown("**Tabla por entrada (EDITABLE)**")

        # Asegura tipos booleanos para que las casillas sean clicables
        for col in ["Flag_Longitud", "Flag_Ambiguedad"]:
            if col in df_entries.columns:
                df_entries[col] = df_entries[col].astype(bool)

        edited_entries = st.data_editor(
            df_entries,
            use_container_width=True,
            height=360,
            num_rows="dynamic",
            column_config={
                "%GC": st.column_config.NumberColumn(
                    "%GC", help="(G+C)/(A+C+G+T) * 100", step=0.001, format="%.3f"
                ),
                "%N": st.column_config.NumberColumn(
                    "%N", help="N/(A+C+G+T+N) * 100", step=0.001, format="%.3f"
                ),
                "Longitud_bases_ACGT": st.column_config.NumberColumn("Longitud ACGT"),
                "Flag_Longitud": st.column_config.CheckboxColumn(
                    "Flag_Longitud", help="Marca si la longitud NO está en el rango 16S."
                ),
                "Flag_Ambiguedad": st.column_config.CheckboxColumn(
                    "Flag_Ambiguedad", help="Marca si %N supera el umbral."
                ),
            },
            key="editor_entries",
        )

        st.download_button(
            "⬇️ Descargar CSV (entradas — con ediciones)",
            edited_entries.to_csv(index=False).encode("utf-8"),
            "gc_entries.csv",
            "text/csv",
        )

        st.markdown("---")
        st.markdown("**Tabla por organismo (resumen)**")
        st.dataframe(df_summary, use_container_width=True, height=320)
        st.download_button(
            "⬇️ Descargar CSV (resumen)",
            df_summary.to_csv(index=False).encode("utf-8"),
            "gc_summary.csv",
            "text/csv",
        )

with tab_grafico:
    st.subheader("3) Gráfico comparativo")
    df_summary = st.session_state.get("df_summary", pd.DataFrame())
    if df_summary.empty:
        st.warning("Sube archivos en la pestaña **Carga**.")
    else:
        opts = list(df_summary["Organismo"])
        subset = st.multiselect("Selecciona organismos a graficar (vacío = todos):", opts, default=opts)
        df_plot = df_summary[df_summary["Organismo"].isin(subset)] if subset else df_summary
        fig = plot_bar(df_plot, title="Contenido de GC en 16S rRNA (promedio por organismo)")
        st.pyplot(fig, use_container_width=True)

        img_bytes = io.BytesIO()
        fig.savefig(img_bytes, format="png", dpi=150, bbox_inches="tight")
        st.download_button("⬇️ Descargar gráfico (PNG)", img_bytes.getvalue(), "gc_barplot.png", "image/png")

with tab_conclusion:
    st.subheader("4) Conclusión automática")
    df_summary = st.session_state.get("df_summary", pd.DataFrame())
    if df_summary.empty:
        st.info("Sube archivos en la pestaña **Carga** para generar la conclusión.")
    else:
        texto = auto_conclusion(df_summary)
        st.markdown(texto)
        st.download_button(
            "⬇️ Descargar conclusión (.txt)",
            texto.encode("utf-8"),
            "conclusion_gc.txt",
            "text/plain"
        )

with tab_guion:
    st.subheader("5) Guion para la presentación (15 min)")

    st.markdown("""
### Acciones que se esperan del agente (según el enunciado)

- Escribir un párrafo que describa los resultados y responda:  
  **¿Cuál organismo tiene mayor cantidad de GC y cuál menor?**  
  → Esto lo realiza automáticamente la pestaña **🧠 Conclusión**.

- Integrar esa funcionalidad en el tablero.  
  → La conclusión se basa en los FASTA cargados en **📥 Carga** y el resumen en **📊 Tabla**.

- Seleccionar cinco microorganismos adicionales de 16S rRNA y agregarlos al conjunto de prueba.  
  → Definidos en `EXTRA_5` y visibles en la descripción de la pestaña **📥 Carga**.

- Plantear fragilidades del código y mitigaciones.  
  → Desarrollado en la pestaña **🧯 Fragilidades**.

- Desarrollar una presentación (retos técnicos, lógica de solución, mejoras).  
  → Guion detallado en esta pestaña **🗂️ Guion**.

- Documentar el código en GitHub.  
  → Este script puede subirse al repositorio con un README que explique entradas, salidas y ejecución.
""")

    guion = f"""
1. Contexto 
   - Objetivo: calcular %GC del 16S rRNA en múltiples bacterias y comparar.
   - Dataset: 9 organismos base + 5 adicionales sugeridos.
   - Entradas: archivos FASTA (16S rRNA); salidas: tabla, gráfico y conclusión automática.

2. Lógica de solución 
   - Lectura de FASTA, limpiando a A/C/G/T/N.
   - Cálculo de %GC sobre ACGT y %N sobre total (ACGTN).
   - Validación de longitud esperada (16S: rango ajustable en la interfaz).
   - Validación de ambigüedad (N): umbral configurable (%N máximo).
   - Agregación por organismo y cálculo de %GC promedio si hay varias entradas.
   - Visualización: tabla editable (por entrada y por organismo) y gráfico de barras.
   - Generación automática del párrafo de conclusión respondiendo quién tiene mayor y menor %GC.

3. Retos técnicos 
   - Encabezados FASTA y nombres de archivo heterogéneos.
   - Secuencias parciales o con alta proporción de N.
   - Varias entradas por organismo → necesidad de agregar y resumir.
   - Rendimiento al manejar muchos archivos y evitar recomputar todo cada vez.
   - Presentar los resultados de forma clara y reutilizable (CSV, PNG, TXT).

4. Resultados clave 
   - Identificación del organismo con mayor %GC y con menor %GC.
   - Rango de variación entre los extremos y promedio general.
   - Interpretación biológica: estabilidad del ADN, genomas reducidos, nichos ecológicos y dependencia del huésped.
   - Relación con bacterias de genoma pequeño y de vida libre.

5. Puntos de mejora 
   - Validación de las secuencias con BLAST 16S para confirmar identidad.
   - Añadir barras de error (±DE) o pruebas estadísticas si hay réplicas.
   - Exportar reporte automático en PDF/HTML desde el mismo dashboard.
   - Crear versión en línea de comandos (CLI) que use los mismos cálculos.
   - Integrar pruebas automáticas (CI) en GitHub para validar cambios de código.

6. Cierre 
   - Demostración rápida: cargar algunos FASTA, mostrar tabla, gráfico y conclusión.
   - Preguntas y discusión sobre posibles extensiones del proyecto.
"""
    st.code(guion.strip())
    st.download_button(
        "⬇️ Descargar guion (.txt)",
        guion.strip().encode("utf-8"),
        "guion_presentacion.txt",
        "text/plain"
    )

with tab_frag:
    st.subheader("6) Fragilidades del código y mitigaciones")
    st.markdown("""
- **Headers / nombres inconsistentes**  
  - Fragilidad: el nombre del organismo se infiere del archivo y puede no coincidir exactamente con NCBI.  
  - Mitigación: permitir editar la tabla y exportar el CSV corregido; documentar la convención de nombres de archivo.

- **Secuencias no 16S o truncadas**  
  - Fragilidad: archivos con regiones parciales o genes diferentes.  
  - Mitigación: usar los umbrales de longitud (min–max) y marcar con `Flag_Longitud` las entradas sospechosas.

- **Bases ambiguas (N)**  
  - Fragilidad: un exceso de N puede sesgar el cálculo de %GC.  
  - Mitigación: calcular %N, usar un umbral configurable y marcar con `Flag_Ambiguedad` para revisar o descartar.

- **Múltiples entradas por organismo**  
  - Fragilidad: resultados inconsistentes si hay secuencias muy distintas.  
  - Mitigación: promediar %GC y %N, y eventualmente calcular desviación estándar o filtrar por calidad.

- **Rendimiento con muchos archivos**  
  - Fragilidad: tiempos de carga largos y consumo de memoria.  
  - Mitigación: uso de caching (`st.cache_data`), procesar por lotes o comprimir entradas (ZIP).

- **Reproducibilidad y despliegue**  
  - Fragilidad: versiones diferentes de librerías o entorno.  
  - Mitigación: incluir `requirements.txt`, especificar versión de Python y documentar comandos de ejecución.

- **Usabilidad / experiencia de usuario**  
  - Fragilidad: que el usuario no comprenda los indicadores de calidad.  
  - Mitigación: tooltips en columnas, mensajes claros en cada pestaña, y guion de presentación que explique la lógica.
""")

