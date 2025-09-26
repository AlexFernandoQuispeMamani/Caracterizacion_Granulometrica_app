# app_streamlit_granulometria.py
# -*- coding: utf-8 -*-
"""
CARACTERIZACIÓN GRANULOMÉTRICA - Web App Streamlit
Autor: Alex Quispe (version corregida)
Instrucciones: instalar requirements.txt y ejecutar:
    streamlit run app_streamlit_granulometria.py
"""
import streamlit as st
import pandas as pd
import numpy as np
import io
import base64
from datetime import datetime
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import qrcode
from PIL import Image
import tempfile
import openpyxl

# ---------------- Configuración ----------------
st.set_page_config(page_title="CARACTERIZACIÓN GRANULOMÉTRICA", layout="wide")

# ---------- Serie Tyler (malla -> abertura en µm) ----------
TYLER = {
    1.05:26500,0.883:22400,0.742:19000,0.624:16000,0.525:13200,0.441:11200,0.371:9500,
    2.5:8000,3:6700,3.5:5600,4:4760,5:4000,6:3350,7:2800,8:2360,9:2000,10:1700,12:1400,
    14:1180,16:1000,20:850,24:710,28:600,32:500,35:425,42:355,48:300,60:250,65:212,
    80:180,100:150,115:125,150:106,170:90,200:75,250:63,270:53,320:45,400:38
}

# ensure unique and sorted mapping (descending by size)
TYLER = {k: v for k, v in sorted(TYLER.items(), key=lambda x: -x[1])}

# ---------- Session state initialization ----------
def ensure_session():
    st.session_state.setdefault('page', 1)
    st.session_state.setdefault('user_info', {})
    st.session_state.setdefault('input_table', pd.DataFrame())
    st.session_state.setdefault('generated_mallas', [])
    st.session_state.setdefault('selected_mode', None)
    st.session_state.setdefault('results_table', pd.DataFrame())
    st.session_state.setdefault('nominal_sizes', pd.DataFrame(columns=['%F(d)','Tamaño (µm)']))
    st.session_state.setdefault('models_fit', {})
    st.session_state.setdefault('uploaded_image', None)
    # We'll store the actual persistent peso total here:
    st.session_state.setdefault('peso_total', 1000.0)
    # helper widget storage
    st.session_state.setdefault('peso_total_input', float(st.session_state.peso_total))

ensure_session()

# ---------- Utilidades ----------
def peso_input_on_change():
    """sincroniza el widget con el valor persistente"""
    try:
        val = float(st.session_state.get('peso_total_input', st.session_state.get('peso_total', 1000.0)))
        st.session_state.peso_total = val
    except Exception:
        pass

def size_label(k):
    return f"{k}# ({TYLER.get(k,'?')} µm)"

def compute_analysis(df_in, mode):
    """Construye la tabla de análisis granulométrico y devuelve DataFrame final (sin formatear)."""
    df = df_in.copy()
    total_weight = float(st.session_state.get('peso_total', 1000.0))

    # Detect column names and normalize
    if 'Abertura (µm)' in df.columns:
        df = df.rename(columns={'Abertura (µm)': 'Tamaño inferior (µm)', 'Peso (g)': 'Peso (g)'})
    if 'Tamaño (µm)' in df.columns and 'Tamaño inferior (µm)' not in df.columns:
        df = df.rename(columns={'Tamaño (µm)': 'Tamaño inferior (µm)'})

    # Ensure necessary columns exist
    if 'Tamaño inferior (µm)' not in df.columns or 'Peso (g)' not in df.columns:
        raise ValueError("La tabla de entrada debe contener columnas 'Tamaño (µm)' y 'Peso (g)' o 'Abertura (µm)' y 'Peso (g)'.")

    df['Tamaño inferior (µm)'] = pd.to_numeric(df['Tamaño inferior (µm)'], errors='coerce')
    df['Peso (g)'] = pd.to_numeric(df['Peso (g)'], errors='coerce').fillna(0.0)

    # ordenar de mayor a menor tamaño inferior (X axis must be decreasing for usual granulometric table)
    df = df.sort_values(by='Tamaño inferior (µm)', ascending=False).reset_index(drop=True)

    # Tamaño superior (para la primera fila se usa multiplicación por sqrt(2))
    size_sup = []
    for i in range(len(df)):
        if i == 0:
            if pd.isna(df.loc[i, 'Tamaño inferior (µm)']):
                sup = np.nan
            else:
                sup = df.loc[i, 'Tamaño inferior (µm)'] * np.sqrt(2)
        else:
            sup = df.loc[i-1, 'Tamaño inferior (µm)']
        size_sup.append(sup)
    df['Tamaño superior (µm)'] = size_sup

    # Calcular peso de la fracción final por diferencia
    suma_pesos = df['Peso (g)'].sum()
    peso_resto = max(total_weight - suma_pesos, 0.0)

    # Añadir última fracción como partículas menores al último tamiz (Tamaño inferior = 0)
    if len(df) > 0:
        last_low = df.loc[len(df)-1, 'Tamaño inferior (µm)']
    else:
        last_low = 0.0
    extra_row = {
        'Tamaño inferior (µm)': 0.0,
        'Tamaño superior (µm)': last_low,
        'Tamaño promedio (µm)': (last_low + 0.0)/2.0,
        'Peso (g)': peso_resto
    }
    df = pd.concat([df, pd.DataFrame([extra_row])], ignore_index=True)

    # Tamaño promedio
    df['Tamaño promedio (µm)'] = (df['Tamaño superior (µm)'] + df['Tamaño inferior (µm)'])/2.0

    # %Peso
    # Evitar división por cero
    if total_weight == 0:
        df['%Peso'] = 0.0
    else:
        df['%Peso'] = 100.0 * df['Peso (g)'] / total_weight

    # %R(d) acumulado retenido y %F(d) pasante
    df['%R(d)'] = df['%Peso'].cumsum()
    df['%F(d)'] = 100.0 - df['%R(d)']

    # Construir fila TOTAL: SOLO mostrar en columnas 'Peso (g)' y '%Peso'
    total_row = {c: np.nan for c in df.columns}
    total_row['Peso (g)'] = total_weight
    # for numeric safety round %Peso sum to 100 (or use df['%Peso'].sum())
    total_row['%Peso'] = round(df['%Peso'].sum(), 6)  # expected to be 100.0
    # If mode SELECT MALLAS, compute Nº de malla intervalos
    cols_order = None
    if mode == "SELECCIONAR MALLAS" and 'Nº Malla (Tyler)' in df_in.columns:
        raw = list(df_in['Nº Malla (Tyler)'])
        intervalos = []
        for i, label in enumerate(raw):
            label = str(label).strip()
            if i == 0:
                intervalos.append(f"{label}")
            else:
                prev = str(raw[i-1]).strip()
                cur = label
                intervalos.append(f"-{prev}+{cur}")
        # append last negative interval -last
        if len(raw) > 0:
            last = str(raw[-1]).strip()
            intervalos.append(f"-{last}")  # última fila se muestra como -NºMalla
        df['Nº de malla (intervalo)'] = intervalos
        total_row['Nº de malla (intervalo)'] = 'Total'

        cols_order = ['Nº de malla (intervalo)', 'Tamaño superior (µm)', 'Tamaño inferior (µm)',
                      'Tamaño promedio (µm)', 'Peso (g)', '%Peso', '%F(d)', '%R(d)']
    else:
        cols_order = ['Tamaño superior (µm)', 'Tamaño inferior (µm)',
                      'Tamaño promedio (µm)', 'Peso (g)', '%Peso', '%F(d)', '%R(d)']

    # Append total row
    df = pd.concat([df, pd.DataFrame([total_row])], ignore_index=True)
    # Reorder columns if possible
    # Some columns may be missing (e.g., 'Nº de malla (intervalo)'); ensure they exist
    for col in cols_order:
        if col not in df.columns:
            df[col] = np.nan
    df = df[cols_order]

    # Format numeric columns to float and round to 2 decimals for display (keep precision in memory)
    numeric_cols = ['Tamaño superior (µm)', 'Tamaño inferior (µm)', 'Tamaño promedio (µm)', 'Peso (g)', '%Peso', '%F(d)', '%R(d)']
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    # return raw numeric DataFrame (rounding applied on display)
    return df

# ---------- PÁGINA 1: Bienvenida ----------
def page_1():
    st.title("CARACTERIZACIÓN GRANULOMÉTRICA")
    st.markdown(f"**Desarrollado por:** Alex Fernando Quispe Mamani — Ingeniero Metalúrgico (UTO).")
    st.markdown("""
    Este programa realiza un análisis granulométrico completo a partir de pesos retenidos por tamiz o
    tamaños ingresados manualmente. Genera tablas, gráficos (varias escalas), estadísticas descriptivas,
    estimaciones según Folk & Ward, ajuste a modelos (GGS, RRSB, Doble Weibull) y exportación de resultados.
    """)
    st.markdown("### INFORMACIÓN GENERAL")
    col1, col2 = st.columns([2,1])

    with col1:
        with st.form("user_info_form", clear_on_submit=False):
            nombre = st.text_input("USUARIO", value=st.session_state['user_info'].get('nombre',''))
            correo = st.text_input("CORREO ELECTRÓNICO", value=st.session_state['user_info'].get('correo',''))
            procedencia = st.text_input("PROCEDENCIA DE LA MUESTRA", value=st.session_state['user_info'].get('procedencia',''))
            codigo = st.text_input("CÓDIGO DE LA MUESTRA", value=st.session_state['user_info'].get('codigo',''))
            fecha_muestreo = st.date_input(
                "FECHA DE MUESTREO",
                value=datetime.fromisoformat(
                    st.session_state['user_info'].get('fecha', datetime.today().date().isoformat())
                ).date()
            )
            inicio = st.form_submit_button("INICIO")

    with col2:
        st.subheader("MET 2260")
        st.image("Imagen2.png", caption="Concentración, Piro y Siderurgia, Adelante Metalurgia!!!...", use_container_width=300)

        if inicio:
            # Validate mandatory fields
            if not nombre or not correo or not procedencia or not codigo:
                st.error("Completa todos los campos obligatorios para continuar.")
            else:
                st.session_state['user_info'] = {
                    'nombre': nombre, 'correo': correo,
                    'procedencia': procedencia, 'codigo': codigo,
                    'fecha': fecha_muestreo.isoformat()
                }
                st.session_state.page = 2
                st.rerun()

# ---------- PÁGINA 2: Datos experimentales ----------
def page_2():
    st.title("DATOS EXPERIMENTALES")
    st.markdown("Inserte tamaños y pesos retenidos sobre cada tamiz para efectuar el análisis granulométrico.")
    st.markdown("Si tiene números de malla use **SELECCIONAR MALLAS**, de lo contrario use **INSERTAR MANUALMENTE**.")

    # Peso total widget (synchronization)
    st.number_input("Peso total (g):", min_value=0.0, value=float(st.session_state.get('peso_total_input', 1000.0)),
                    step=0.1, key='peso_total_input', on_change=peso_input_on_change)
    # garantizamos que la variable persistente esté sincronizada
    st.session_state.peso_total = float(st.session_state.get('peso_total', st.session_state.get('peso_total_input', 1000.0)))

    st.markdown("**Modo de inserción de datos**")
    mode = st.radio("Selecciona modo:", ["SELECCIONAR MALLAS", "INSERTAR MANUALMENTE"], index=0, key='mode_radio')
    st.session_state.selected_mode = mode

    if mode == "SELECCIONAR MALLAS":
        st.info("Selecciona las mallas de la serie Tyler. Luego pulsa **Generar tabla de mallas**.")
        # build labels for multiselect
        malla_items = sorted(TYLER.items(), key=lambda x: -x[1])
        labels = [f"{int(k) if float(k).is_integer() else k}# - {v} µm" for k, v in malla_items]
        selected_labels = st.multiselect("Selecciona mallas (múltiple, ordenadas descendente si quieres):", labels, default=st.session_state.get('generated_mallas_labels', []))

        if st.button("Generar tabla de mallas"):
            # convert labels to keys
            selected_keys = []
            rows = []
            for lab in selected_labels:
                key_str = lab.split('#')[0]
                try:
                    k = int(key_str)
                except:
                    try:
                        k = float(key_str)
                    except:
                        k = key_str
                selected_keys.append(k)
            st.session_state['generated_mallas'] = selected_keys
            st.session_state['generated_mallas_labels'] = selected_labels

            for k in selected_keys:
                rows.append({'Nº Malla (Tyler)': str(k)+'#', 'Abertura (µm)': TYLER.get(k, np.nan), 'Peso (g)': np.nan})
            if len(rows) == 0:
                st.warning("Selecciona al menos una malla.")
            else:
                st.session_state.input_table = pd.DataFrame(rows)
                st.success("Tabla generada. Completa la columna 'Peso (g)' con tus datos y pulsa EJECUTAR.")
                st.rerun()

    else:  # INSERTAR MANUALMENTE
        st.info("Insertar manualmente los tamaños (µm) y pesos (g).")
        n = st.number_input("Número de filas a insertar (3-25):", min_value=3, max_value=25, value=6, step=1, key='n_rows_input')
        if st.button("Generar tabla manual"):
            df = pd.DataFrame({'Tamaño (µm)': [np.nan]*int(n), 'Peso (g)': [np.nan]*int(n)})
            st.session_state.input_table = df
            st.success("Tabla generada. Completa los tamaños y pesos y pulsa EJECUTAR.")
            st.rerun()

    st.markdown("**Tabla de entrada** (edítala y luego pulsa EJECUTAR):")
    if not st.session_state.input_table.empty:
        edited = st.data_editor(st.session_state.input_table, num_rows="dynamic", key='input_table_editor')
        st.session_state.input_table = edited

    col1, col2, col3 = st.columns([1,1,1])
    with col1:
        if st.button("ANTERIOR"):
            st.session_state.page = 1
            st.rerun()
    with col2:
        if st.button("EJECUTAR"):
            # validation
            if st.session_state.input_table.empty:
                st.error("Genera y completa la tabla antes de ejecutar.")
            else:
                # compute and store results_table ONCE (this prevents recompute on widget change)
                try:
                    results = compute_analysis(st.session_state.input_table, st.session_state.selected_mode)
                    st.session_state.results_table = results
                    st.success("Análisis calculado correctamente.")
                    st.session_state.page = 3
                    st.rerun()
                except Exception as e:
                    st.error(f"Error al calcular: {e}")
    with col3:
        st.write("")  # spacing

# ---------- PÁGINA 3: Análisis granulométrico (Resultados) ----------
def page_3():
    st.title("ANÁLISIS GRANULOMÉTRICO")
    results = st.session_state.results_table.copy()
    if results.empty:
        st.error("No hay resultados calculados. Vuelve a la página de Datos Experimentales y pulsa EJECUTAR.")
        if st.button("Regresar a DATOS EXPERIMENTALES"):
            st.session_state.page = 2
            st.rerun()
        return

    # Mostrar tabla de resultados
    fmt = {}
    for c in results.columns:
        if c in ['Peso (g)', '%Peso', '%F(d)', '%R(d)', 'Tamaño promedio (µm)', 'Tamaño inferior (µm)', 'Tamaño superior (µm)']:
            fmt[c] = "{:.2f}"
        else:
            fmt[c] = "{}"
    st.markdown("**Tabla de resultados** (última fila = TOTAL).")
    st.dataframe(results.style.format(fmt), height=320)

    # Selección de gráfico y escala
    st.markdown("**Seleccione gráfico**")
    grafico = st.selectbox("SELECCIONE GRÁFICO", [
        "Histograma de frecuencia",
        "Diagrama de simple distribución",
        "Diagrama Acumulativo de Subtamaño",
        "Diagrama Acumulativo de Sobretamaño",
        "Diagrama Acumulativo (Combinación)",
        "Curvas granulométricas (Combinación 2,3,4)"
    ], index=0, key='grafico_select')

    escala = st.selectbox("Escala", ["Escala decimal", "Escala semilogarítmica (X log)", "Escala logarítmica (ambos log)"], index=0, key='escala_select')

    # Usar Tamaño inferior como X en curvas acumulativas y de distribución
    plot_df = results.iloc[:-1].copy()  # quitar fila TOTAL
    plot_df = plot_df[plot_df['Tamaño inferior (µm)'].notna()]

    x_inf = plot_df['Tamaño inferior (µm)'].replace(0, np.nan)
    x_avg = plot_df['Tamaño promedio (µm)'].replace(0, np.nan)  # solo para histograma
    y_pct = plot_df['%Peso']
    yf = plot_df['%F(d)']
    yr = plot_df['%R(d)']

    fig, ax = plt.subplots(figsize=(9, 4))
    fig.patch.set_facecolor('#e6e6e6')
    ax.set_facecolor('white')

    lw = 0.9
    ms = 4

    if grafico == "Histograma de frecuencia":
        ax.bar(
            x_avg, y_pct,
            width=(np.nanmax(x_avg)/len(x_avg) if len(x_avg) > 0 else 1),
            color="black", edgecolor="black", linewidth=0.5, alpha=0.9
        )
        ax.set_xlabel("Tamaño promedio (µm)")
        ax.set_ylabel("%Peso")

    elif grafico == "Diagrama de simple distribución":
        ax.plot(x_inf, y_pct, marker='o', markersize=ms, linewidth=lw, color='black')
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%Peso")

    elif grafico == "Diagrama Acumulativo de Subtamaño":
        ax.plot(x_inf, yf, marker='o', markersize=ms, linewidth=lw, color='black')
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%F(d)")

    elif grafico == "Diagrama Acumulativo de Sobretamaño":
        ax.plot(x_inf, yr, marker='x', markersize=ms, linewidth=lw, color='black')
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%R(d)")

    elif grafico == "Diagrama Acumulativo (Combinación)":
        ax.plot(x_inf, yf, label='%F(d)', marker='o', markersize=ms, linewidth=lw, color='black')
        ax.plot(x_inf, yr, label='%R(d)', marker='x', markersize=ms, linewidth=lw, color='black')
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("Porcentaje")
        ax.legend()

    else:  # Curvas granulométricas (Combinación 2,3,4)
        ax.plot(x_inf, y_pct, label='%Peso', marker='s', markersize=ms, linewidth=lw, color='black')
        ax.plot(x_inf, yf, label='%F(d)', marker='o', markersize=ms, linewidth=lw, color='black')
        ax.plot(x_inf, yr, label='%R(d)', marker='x', markersize=ms, linewidth=lw, color='black')
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("Porcentaje")
        ax.legend()

    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 10))

    # Escalas
    if escala == "Escala semilogarítmica (X log)":
        ax.set_xscale('log')
    elif escala == "Escala logarítmica (ambos log)":
        ax.set_xscale('log')
        ax.set_yscale('log')

    ax.grid(True, which='both', ls='--', alpha=0.5)
    ax.set_title(grafico)
    st.pyplot(fig)

    # Navegación
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            st.session_state.page = 2
            st.rerun()
    with col2:
        if st.button("SIGUIENTE"):
            st.session_state.page = 4
            st.rerun()

# ---------- PÁGINA 4: Análisis de datos ----------
def page_4():
    st.title("ANÁLISIS DE DATOS")
    results = st.session_state.results_table.copy()
    if results.empty:
        st.error("No hay resultados calculados. Vuelve a DATOS EXPERIMENTALES y pulsa EJECUTAR.")
        return

    # Sección 1: ANÁLISIS ESTADÍSTICO
    st.header("ANÁLISIS ESTADÍSTICO")
    sizes = results['Tamaño inferior (µm)'].iloc[:-1].replace(0, np.nan).dropna()
    weights = results['%Peso'].iloc[:-1].replace(0, np.nan).dropna() / 100.0
    mask = (~sizes.isna()) & (~weights.isna())
    sizes = sizes.loc[mask]
    weights = weights.loc[mask]

    if len(sizes) == 0:
        st.warning("No hay suficientes datos para el análisis estadístico.")
    else:
        wsum = weights.sum()
        if wsum == 0:
            mean = float(sizes.mean())
        else:
            mean = float((sizes * weights).sum() / wsum)
        try:
            median = float(np.interp(50, 100*weights.cumsum(), sizes))
        except:
            median = float(np.median(sizes))
        mode_idx = np.argmax(results['%Peso'].iloc[:-1].values)
        mode = float(results['Tamaño inferior (µm)'].iloc[:-1].values[mode_idx])
        variance = float(((sizes - mean)**2 * weights).sum() / (weights.sum() if weights.sum()>0 else len(sizes)))
        std = float(np.sqrt(variance))
        kurtosis = float(pd.Series(sizes).kurtosis())
        rango = float(np.nanmax(sizes) - np.nanmin(sizes))

        stats_tbl = pd.DataFrame({
            'Estadístico':['Media (µm)','Mediana (µm)','Moda (µm)','Varianza (µm²)',
                           'Desvío estándar (µm)','Curtosis','Rango (µm)'],
            'Valor':[mean, median, mode, variance, std, kurtosis, rango]
        })
        colA, colB = st.columns([1,1])
        with colA:
            st.table(stats_tbl.style.format({ 'Valor':'{:.3f}'}))
        with colB:
            comments = []
            if variance < (mean * 0.01 if mean>0 else 1):
                comments.append(f"Varianza {variance:.3f}: variabilidad de tamaños **muy baja**; partículas relativamente uniformes.")
            elif variance < (mean * 0.1 if mean>0 else 10):
                comments.append(f"Varianza {variance:.3f}: variabilidad **moderada** en tamaños.")
            else:
                comments.append(f"Varianza {variance:.3f}: variabilidad **alta** en tamaños (amplia dispersión granulométrica).")

            if rango <= 50:
                comments.append(f"Rango {rango:.1f} µm: espectro granulométrico **estrecho** (partículas similares).")
            elif rango <= 300:
                comments.append(f"Rango {rango:.1f} µm: espectro **moderado**.")
            else:
                comments.append(f"Rango {rango:.1f} µm: espectro **amplio** (muestra heterogénea).")

            if kurtosis > 1.5:
                comments.append(f"Curtosis {kurtosis:.2f}: distribución **leptocúrtica** (cola delgada, pico alto).")
            elif kurtosis < -1.0:
                comments.append(f"Curtosis {kurtosis:.2f}: distribución **platocúrtica** (colas gruesas, pico bajo).")
            else:
                comments.append(f"Curtosis {kurtosis:.2f}: distribución **mesocúrtica** (casi normal).")

            st.markdown("**Interpretación automática:**")
            for c in comments:
                st.write("- " + c)

    # Sección 2: Registro de tamaños y porcentajes pasantes
    st.header("REGISTRO DE TAMAÑOS Y PORCENTAJES PASANTES")
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("PERFIL GRANULÓMETRICO")
        st.markdown("El perfil granulométrico de una muestra describe cualitativamente la distribución por tamaños de dicha muestra.")
        plot_df = results.iloc[:-1].copy()
        plot_df = plot_df[plot_df['Tamaño inferior (µm)'].notna()]
        x = plot_df['Tamaño inferior (µm)']
        y = plot_df['%F(d)']

        fig, ax = plt.subplots(figsize=(6,4))
        fig.patch.set_facecolor('#e6e6e6')
        ax.set_facecolor('white')
        mask_plot = x.notna() & y.notna()
        x_plot = x[mask_plot]
        y_plot = y[mask_plot]

        ax.plot(x_plot, y_plot, color='black', linewidth=0.8, zorder=2)
        ax.scatter(x_plot, y_plot, facecolors='white', edgecolors='black', s=30,
                   linewidths=0.8, zorder=3, marker='o')

        ax.set_xlabel("Tamaño (µm)")
        ax.set_ylabel("%F(d)")
        ax.set_ylim(0, 100)
        ax.set_yticks(np.arange(0, 101, 10))
        ax.grid(True, which='both', ls='--', alpha=0.5)
        st.pyplot(fig)

        try:
            mask = (~np.isnan(x)) & (~np.isnan(y))
            xs = np.array(x[mask])
            ys = np.array(y[mask])
            order = np.argsort(xs)
            xs_s = xs[order]
            ys_s = ys[order]
            inv = interp1d(ys_s, xs_s, fill_value="extrapolate", bounds_error=False)
            d80 = float(inv(80.0))
            st.markdown(f"**d80 = {d80:.2f} µm** — El 80% de partículas es menor a este tamaño.")
        except Exception as e:
            st.warning("No hay suficientes datos para calcular d80: " + str(e))

    with col2:
        st.subheader("%F(d) y TAMAÑOS NOMINALES")
        modo_calc = st.selectbox("Seleccione operación:", ["CALCULAR TAMAÑO", "CALCULAR %F(d)"], key='modo_calc')
        if modo_calc == "CALCULAR TAMAÑO":
            pct_input = st.number_input("Ingrese %F(d) =", min_value=0.0, max_value=100.0,
                                        value=50.0, key='pct_input')
            if st.button("GRABAR d (pct->tamaño)"):
                x = plot_df['Tamaño inferior (µm)']
                y = plot_df['%F(d)']
                mask = (~np.isnan(x)) & (~np.isnan(y))
                try:
                    inv = interp1d(y[mask], x[mask], fill_value="extrapolate", bounds_error=False)
                    d_val = float(inv(pct_input))
                    st.session_state.nominal_sizes = pd.concat([
                        st.session_state.nominal_sizes,
                        pd.DataFrame([{'%F(d)': pct_input, 'Tamaño (µm)': d_val}])
                    ], ignore_index=True)
                    st.success(f"Grabado: {pct_input}% -> {d_val:.2f} µm")
                except Exception as e:
                    st.error("No se puede interpolar: " + str(e))
        else:
            d_input = st.number_input("Ingrese d (µm) =", min_value=0.0, value=100.0, key='d_input')
            if st.button("GRABAR %F(d) (tamaño->pct)"):
                x = plot_df['Tamaño inferior (µm)']
                y = plot_df['%F(d)']
                mask = (~np.isnan(x)) & (~np.isnan(y))
                try:
                    f = interp1d(x[mask], y[mask], fill_value="extrapolate", bounds_error=False)
                    pctv = float(f(d_input))
                    st.session_state.nominal_sizes = pd.concat([
                        st.session_state.nominal_sizes,
                        pd.DataFrame([{'%F(d)': pctv, 'Tamaño (µm)': d_input}])
                    ], ignore_index=True)
                    st.success(f"Grabado: {d_input} µm -> {pctv:.2f}%")
                except Exception as e:
                    st.error("No se puede interpolar: " + str(e))

        st.markdown("**Tabla de tamaños nominales grabados (editable)**")
        if not st.session_state.nominal_sizes.empty:
            st.session_state.nominal_sizes = st.data_editor(
                st.session_state.nominal_sizes,
                num_rows="dynamic"
            )
        else:
            st.write("_Aún no hay tamaños grabados._")

    # Sección 3: Folk & Ward
    st.header("ESTADÍSTICOS SEGÚN FOLK & WARD")
    st.markdown("""
    **Folk & Ward (1957):** se calculan parámetros de tendencia central y dispersión (M, Md, σ, Sk, K).  
    Se utiliza la **escala phi (φ)**, que corresponde a una transformación logarítmica negativa de los tamaños de partícula en milímetros:

    φ = –log₂(d [mm])  

    Esta escala permite comparar distribuciones granulométricas de manera estandarizada, ya que convierte los tamaños (que cubren varios órdenes de magnitud) en una escala lineal en términos estadísticos.  
    """)

    if st.button("ESTIMAR"):
        try:
            req_pcts = [5,16,25,50,75,84,95]

            if st.session_state.nominal_sizes.empty:
                st.warning("Grabe los tamaños nominales d5, d16, d25, d50, d75, d84, d95")
            else:
                ds = {}
                for p in req_pcts:
                    match = st.session_state.nominal_sizes.loc[
                        (st.session_state.nominal_sizes['%F(d)'].sub(p).abs() < 1e-6), 'Tamaño (µm)']
                    if match.empty:
                        st.warning(f"Falta grabar d{p}")
                        ds = None
                        break
                    else:
                        ds[p] = float(match.iloc[0])

                if ds is not None:
                    # Reordenar por percentil para evitar problemas de grabación desordenada
                    ds = dict(sorted(ds.items(), key=lambda kv: kv[0]))

                    # Convertir a mm y luego a escala phi
                    ds_mm = {p: ds[p]/1000.0 for p in ds}
                    for p in ds_mm:
                        if ds_mm[p] <= 0:
                            raise ValueError("Valores de tamaño no positivos para convertir a phi.")
                    phi = {p: -np.log2(ds_mm[p]) for p in ds_mm}

                    # Cálculos Folk & Ward (usando fórmulas manuales)
                    M = (phi[16] + phi[50] + phi[84]) / 3.0
                    Md = phi[50]
                    sigmaI = ((phi[84] - phi[16]) / 4.0) + ((phi[95] - phi[5]) / 6.6)
                    SkI = ((phi[16] + phi[84] - 2*phi[50]) / (2*(phi[84]-phi[16]))) + \
                          ((phi[5] + phi[95] - 2*phi[50]) / (2*(phi[95]-phi[5])))
                    KG = (phi[95] - phi[5]) / (2.44 * (phi[75] - phi[25]))

                    # Tabla sin índice
                    folk_tbl = pd.DataFrame({
                        'Parámetro':['M (φ)','Md (φ)','σ (φ)','Sk (φ)','K (φ)'],
                        'Valor':[M, Md, sigmaI, SkI, KG]
                    })

                    col_tbl, col_interp = st.columns([1,1])
                    with col_tbl:
                        st.subheader("Tabla Folk & Ward")
                        mostrar_indice = st.checkbox("Mostrar índice", value=True)
                        st.table(folk_tbl.style.format({'Valor':'{:.3f}'}), hide_index=not mostrar_indice)

                    with col_interp:
                        st.subheader("Interpretación Folk & Ward")

                        # Media
                        if M < 0:
                            mean_comment = f"M = {M:.3f} φ → predominan partículas muy gruesas."
                        elif M < 1:
                            mean_comment = f"M = {M:.3f} φ → predominan partículas gruesas."
                        elif M < 2:
                            mean_comment = f"M = {M:.3f} φ → predominan partículas de tamaño medio."
                        elif M < 3:
                            mean_comment = f"M = {M:.3f} φ → predominan partículas finas."
                        elif M < 4:
                            mean_comment = f"M = {M:.3f} φ → predominan partículas muy finas."
                        else:
                            mean_comment = f"M = {M:.3f} φ → predominan finos (limo o arcilla)."

                        # Dispersión (σ)
                        if sigmaI < 0.35:
                            sigma_comment = f"σ = {sigmaI:.3f} → granulometría homogénea, distribución uniforme (espectro muy estrecho)."
                        elif sigmaI < 0.50:
                            sigma_comment = f"σ = {sigmaI:.3f} → granulometría relativamente homogénea (espectro estrecho)."
                        elif sigmaI < 0.71:
                            sigma_comment = f"σ = {sigmaI:.3f} → granulometría moderada (espectro moderadamente estrecho)."
                        elif sigmaI < 1.00:
                            sigma_comment = f"σ = {sigmaI:.3f} → granulometría heterogénea moderada (espectro intermedio)."
                        elif sigmaI < 2.00:
                            sigma_comment = f"σ = {sigmaI:.3f} → granulometría heterogénea (espectro amplio)."
                        else:
                            sigma_comment = f"σ = {sigmaI:.3f} → granulometría muy heterogénea (espectro muy amplio)."

                        # Sesgo
                        if SkI > 0:
                            skew_comment = f"Sk = {SkI:.3f} → distribución sesgada hacia finos (exceso de partículas pequeñas)."
                        elif SkI < 0:
                            skew_comment = f"Sk = {SkI:.3f} → distribución sesgada hacia gruesos (exceso de partículas grandes)."
                        else:
                            skew_comment = f"Sk = {SkI:.3f} → distribución aproximadamente simétrica."

                        # Curtosis
                        if KG < 0.67:
                            kurt_comment = f"K = {KG:.3f} → distribución platicúrtica (aplanada, mezcla amplia de tamaños)."
                        elif KG < 1.11:
                            kurt_comment = f"K = {KG:.3f} → distribución mesocúrtica (curva normal, típica de molienda estándar)."
                        else:
                            kurt_comment = f"K = {KG:.3f} → distribución leptocúrtica (pico agudo, concentración de tamaños, control granulométrico preciso)."

                        # Mostrar interpretaciones
                        st.write("- " + mean_comment)
                        st.write("- " + sigma_comment)
                        st.write("- " + skew_comment)
                        st.write("- " + kurt_comment)

        except Exception as e:
            st.warning("No se pueden calcular Folk & Ward: " + str(e))

    # Navegación
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            st.session_state.page = 3
            st.rerun()
    with col2:
        if st.button("SIGUIENTE"):
            st.session_state.page = 5
            st.rerun()

            
# ---------- MODELOS (GGS, RRSB, Doble Weibull) ----------
def GGS_model(d, m, Dm):
    # flexible cumulative form: 100/(1+(d/Dm)**(-m))
    # avoid division by zero
    d = np.array(d, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return 100.0 * (1.0 / (1.0 + (d/Dm)**(-m)))

def RRSB_model(d, m, l):
    d = np.array(d, dtype=float)
    return 100.0 * (1 - np.exp(- (d/l)**m))

def double_weibull(d, alpha, k1, k2, d80):
    # d1,d2 heuristic from d80
    d1 = d80 * 0.6
    d2 = d80 * 1.4
    d = np.array(d, dtype=float)
    return 100.0 * (alpha*(1 - np.exp(-(d/d1)**k1)) + (1-alpha)*(1 - np.exp(-(d/d2)**k2)))

# ---------- PÁGINA 5: Selección del Modelo ----------
def page_5():
    st.title("SELECCIÓN DEL MODELO")
    st.markdown("Ajuste de modelos: GGS, RRSB y Doble Weibull. Se estiman parámetros minimizando SSE (F.O.).")
    results = st.session_state.results_table.copy()
    if results.empty:
        st.error("No hay resultados para ajustar. Regresa a DATOS EXPERIMENTALES y pulsa EJECUTAR.")
        return

    # Prepare data for fitting (excluir fila TOTAL)
    df_fit = results.iloc[:-1].copy()
    mask = (df_fit['Tamaño inferior (µm)'] > 0) & (~np.isnan(df_fit['%F(d)']))
    d = df_fit['Tamaño inferior (µm)'][mask].astype(float).values   # usar tamaño inferior
    y_exp = df_fit['%F(d)'][mask].astype(float).values

    if len(d) < 3:
        st.warning("Se requieren al menos 3 puntos válidos para ajuste.")
        st.write("Si necesitas más puntos, vuelve a DATOS EXPERIMENTALES y completa la tabla.")
        return

    if st.button("AJUSTAR"):
        # GGS
        try:
            x0 = [1.0, np.median(d)]
            def f_ggs(params):
                m, Dm = params
                ypred = GGS_model(d, m, Dm)
                return np.sum((y_exp - ypred)**2)
            res1 = minimize(f_ggs, x0, bounds=[(0.01,10),(1e-6, max(d)*10)])
            FO_ggs = float(res1.fun)
            ggs_params = res1.x.tolist()
        except Exception:
            FO_ggs = np.inf; ggs_params = [np.nan, np.nan]

        # RRSB
        try:
            x0 = [1.0, np.median(d)]
            def f_rrsb(params):
                m, l = params
                ypred = RRSB_model(d, m, l)
                return np.sum((y_exp - ypred)**2)
            res2 = minimize(f_rrsb, x0, bounds=[(0.01,10),(1e-6,max(d)*10)])
            FO_rrsb = float(res2.fun)
            rrsb_params = res2.x.tolist()
        except Exception:
            FO_rrsb = np.inf; rrsb_params = [np.nan, np.nan]

        # Double Weibull
        try:
            try:
                inv = interp1d(df_fit['%F(d)'], df_fit['Tamaño inferior (µm)'], 
                               fill_value="extrapolate", bounds_error=False)
                init_d80 = float(inv(80.0))
                if init_d80 <= 0 or np.isnan(init_d80):
                    init_d80 = np.median(d)
            except:
                init_d80 = np.median(d)
            x0 = [0.5, 1.0, 1.0, init_d80]
            bounds_dw = [(0.0,1.0),(0.01,10.0),(0.01,10.0),(1e-3,max(d)*10)]
            def f_double(params):
                alpha, k1, k2, d80 = params
                ypred = double_weibull(d, alpha, k1, k2, d80)
                return np.sum((y_exp - ypred)**2)
            res3 = minimize(f_double, x0, bounds=bounds_dw)
            FO_dw = float(res3.fun)
            dw_params = res3.x.tolist()
        except Exception:
            FO_dw = np.inf; dw_params = [np.nan, np.nan, np.nan, np.nan]

        # Save fits
        st.session_state.models_fit = {
            'GGS': {'FO':FO_ggs, 'params':ggs_params},
            'RRSB': {'FO':FO_rrsb, 'params':rrsb_params},
            'DoubleWeibull': {'FO':FO_dw, 'params':dw_params}
        }
        st.success("Ajustes completados.")

    # ----------- Mostrar tabla comparativa ----------
    if st.session_state.get("models_fit", None):
        fits = st.session_state.models_fit

        table_data = []
        d_inferior = df_fit['Tamaño inferior (µm)'][mask].astype(float).values

        for i in range(len(d_inferior)):
            row = {
                "Tamaño inferior (µm)": d_inferior[i],
                "%F(d)e": y_exp[i]
            }
            # GGS
            if not np.isnan(fits['GGS']['params']).any():
                m, Dm = fits['GGS']['params']
                y_g = GGS_model([d[i]], m, Dm)[0]
                row["%F(d)m_GGS"] = y_g
                row["ε²_GGS"] = (y_exp[i]-y_g)**2
            else:
                row["%F(d)m_GGS"] = np.nan
                row["ε²_GGS"] = np.nan

            # RRSB
            if not np.isnan(fits['RRSB']['params']).any():
                m2, l = fits['RRSB']['params']
                y_r = RRSB_model([d[i]], m2, l)[0]
                row["%F(d)m_RRSB"] = y_r
                row["ε²_RRSB"] = (y_exp[i]-y_r)**2
            else:
                row["%F(d)m_RRSB"] = np.nan
                row["ε²_RRSB"] = np.nan

            # Double Weibull
            if not np.isnan(fits['DoubleWeibull']['params']).any():
                alpha, k1, k2, d80 = fits['DoubleWeibull']['params']
                y_dw = double_weibull([d[i]], alpha, k1, k2, d80)[0]
                row["%F(d)m_DW"] = y_dw
                row["ε²_DW"] = (y_exp[i]-y_dw)**2
            else:
                row["%F(d)m_DW"] = np.nan
                row["ε²_DW"] = np.nan

            table_data.append(row)

        df_comp = pd.DataFrame(table_data)
        df_comp = df_comp[[
            "Tamaño inferior (µm)", "%F(d)e",
            "%F(d)m_GGS","ε²_GGS",
            "%F(d)m_RRSB","ε²_RRSB",
            "%F(d)m_DW","ε²_DW"
        ]]

        st.subheader("Tabla comparativa: Experimental vs Modelos")
        st.dataframe(df_comp.style.format("{:.2f}"), height=320)

    # ----------- Comparación de FO y gráficas -----------
    if st.session_state.models_fit:
        fits = st.session_state.models_fit
        fo_tbl = pd.DataFrame([
            {'Modelo':'GGS','F.O.':fits['GGS']['FO']},
            {'Modelo':'RRSB','F.O.':fits['RRSB']['FO']},
            {'Modelo':'Doble Weibull','F.O.':fits['DoubleWeibull']['FO']}
        ])
        st.subheader("Comparación de funciones objetivo (F.O.)")
        st.table(fo_tbl.style.format({'F.O.':'{:.6g}'}))

        best = min(fits.items(), key=lambda x: x[1]['FO'])
        best_model_name = best[0]
        st.markdown(f"**Mejor modelo:** {best_model_name} con F.O. = {best[1]['FO']:.6g}")

        # Graficar
        df_fit_plot = results.iloc[:-1].copy()
        df_fit_plot = df_fit_plot[df_fit_plot['Tamaño inferior (µm)'].notna()]
        xdata = df_fit_plot['Tamaño inferior (µm)'].values
        ydata = df_fit_plot['%F(d)'].values

        fig, ax = plt.subplots(figsize=(8,4))
        ax.plot(xdata, ydata, 'o', label='Experimental', color='black', markersize=4)
        dd = np.linspace(np.min(xdata), np.max(xdata), 300)
        # GGS
        if not np.isnan(fits['GGS']['params']).any():
            m, Dm = fits['GGS']['params']
            ax.plot(dd, GGS_model(dd, m, Dm), '-', label=f'GGS (m={m:.3f}, Dm={Dm:.3f})', linewidth=0.9)
        # RRSB
        if not np.isnan(fits['RRSB']['params']).any():
            m2, l = fits['RRSB']['params']
            ax.plot(dd, RRSB_model(dd, m2, l), '--', label=f'RRSB (m={m2:.3f}, l={l:.3f})', linewidth=0.9)
        # Double Weibull
        p = fits['DoubleWeibull']['params']
        if not np.isnan(p).any():
            alpha, k1, k2, d80 = p
            ax.plot(dd, double_weibull(dd, alpha, k1, k2, d80), ':', label=f'DW (α={alpha:.3f},k1={k1:.3f},k2={k2:.3f},d80={d80:.2f})', linewidth=0.9)
        ax.set_xlabel("Tamaño (µm)")
        ax.set_ylabel("%F(d)")
        ax.set_ylim(0, 100)
        ax.grid(True, ls='--', alpha=0.5)
        ax.legend()
        st.pyplot(fig)

        if best[1]['FO'] > 1e6:
            st.warning("Ningún modelo representa bien los datos experimentales (F.O. muy grande).")
        else:
            st.success(f"El mejor modelo que representa los resultados experimentales es {best_model_name} (F.O. = {best[1]['FO']:.6g}).")

    # Navegación
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            st.session_state.page = 4
            st.rerun()  
    with col2:
        if st.button("SIGUIENTE"):
            st.session_state.page = 6
            st.rerun()

# ---------- PÁGINA 6: Exportación ----------
def page_6():
    st.title("EXPORTACIÓN DE DATOS")
    st.markdown("Descargar todas las tablas en un archivo Excel. Cada usuario obtiene su propio archivo en modo solo lectura.")

    # Construir Excel en memoria
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # Info usuario
        if "user_info" in st.session_state:
            ui = pd.DataFrame([st.session_state['user_info']])
            ui.to_excel(writer, sheet_name='InfoUsuario', index=False)
        # Tabla de entrada
        if "input_table" in st.session_state and not st.session_state.input_table.empty:
            st.session_state.input_table.to_excel(writer, sheet_name='Entrada', index=False)
        # Resultados
        if "results_table" in st.session_state and not st.session_state.results_table.empty:
            st.session_state.results_table.to_excel(writer, sheet_name='Resultados', index=False)
        # Tamaños nominales
        if "nominal_sizes" in st.session_state and not st.session_state.nominal_sizes.empty:
            st.session_state.nominal_sizes.to_excel(writer, sheet_name='TamañosNominales', index=False)
        # Modelos
        if "models_fit" in st.session_state and st.session_state.models_fit:
            models_df = pd.DataFrame([
                {'Modelo': k, 'FO': v['FO'], 'Parametros': str(v['params'])}
                for k, v in st.session_state.models_fit.items()
            ])
            models_df.to_excel(writer, sheet_name='Modelos', index=False)

    data = output.getvalue()

    # Botón oficial de Streamlit para descargar
    st.download_button(
        label="📥 Descargar Excel",
        data=data,
        file_name="analisis_granulometrico.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    )

    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            st.session_state.page = 5
            st.rerun()
    with col2:
        if st.button("VOLVER AL INICIO"):
            st.session_state.page = 1
            st.rerun()

# ---------- Router ----------
def main():
    ensure_session()
    page = st.session_state.page
    if page == 1:
        page_1()
    elif page == 2:
        page_2()
    elif page == 3:
        page_3()
    elif page == 4:
        page_4()
    elif page == 5:
        page_5()
    elif page == 6:
        page_6()

if __name__ == '__main__':
    main()



























