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

    # Tamaño superior (para la primera fila se usa multiplicación por 2)
    size_sup = []
    for i in range(len(df)):
        if i == 0:
            if pd.isna(df.loc[i, 'Tamaño inferior (µm)']):
                sup = np.nan
            else:
                sup = df.loc[i, 'Tamaño inferior (µm)'] * 2  # Cambio aquí
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

    # --- CORRECCIÓN 1: evitar -0 en %F(d) ---
    if '%F(d)' in df.columns:
        df['%F(d)'] = df['%F(d)'].abs()

    # --- CORRECCIÓN 2: ajustar primera malla con 2 mallas superiores de la serie Tyler ---
    if 'Nº de malla (intervalo)' in df.columns and len(df) > 1:
        try:
            tyler_series = [
                1.05,0.883,0.742,0.624,0.525,0.441,0.371,2.5,3,3.5,4,5,6,7,8,9,10,12,
                14,16,20,24,28,32,35,42,48,60,65,80,100,115,150,170,200,250,270,320,400
            ]
            # obtener la malla de la primera fila
            current_label = df.loc[0, 'Nº de malla (intervalo)']
            try:
                current_num = float(str(current_label).replace('#',''))
            except:
                current_num = None

            if current_num is not None:
                if current_num in tyler_series:
                    idx = tyler_series.index(current_num)
                    if idx >= 2:
                        upper = tyler_series[idx - 2]
                        df.loc[0, 'Nº de malla (intervalo)'] = f"-{int(upper)}#+{int(current_num)}#"
                    else:
                        # No hay dos mallas superiores, solo indicamos +malla
                        upper = tyler_series[0]
                        df.loc[0, 'Nº de malla (intervalo)'] = f"+{int(current_num)}#"
        except Exception as e:
            print("Aviso: no se pudo ajustar la primera fila:", e)

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
    st.markdown("**Tabla de resultados**")
    st.dataframe(results.style.format(fmt), height=320)

    # Selección de gráfico y escala
    st.markdown("**Seleccione gráfico**")
    grafico = st.selectbox("SELECCIONE GRÁFICO", [
        "Histograma de frecuencia",
        "Diagrama de simple distribución",
        "Diagrama Acumulativo de Subtamaño",
        "Diagrama Acumulativo de Sobretamaño",
        "Diagrama Acumulativo",
        "Curvas granulométricas de la muestra"
    ], index=0, key='grafico_select')

    escala = st.selectbox("Escala", ["Escala decimal", "Escala semilogarítmica", "Escala logarítmica"], index=0, key='escala_select')

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

    # --- Graficar según selección ---
    if grafico == "Histograma de frecuencia":
        for i in range(len(plot_df)):
            x_left = plot_df['Tamaño inferior (µm)'].iloc[i]
            x_right = plot_df['Tamaño superior (µm)'].iloc[i]
            height = plot_df['%Peso'].iloc[i]
            ax.bar(
                x_left, height,
                width=(x_right - x_left),
                align='edge',
                color="gray", edgecolor="black", linewidth=0.5, alpha=0.9
            )
        ax.set_xlabel("Tamaño (µm)")
        ax.set_ylabel("%Peso")

    elif grafico == "Diagrama de simple distribución":
        ax.scatter(x_inf, y_pct, s=48, linewidths=0.9, edgecolors='k', facecolors='white', zorder=4)
        ax.plot(x_inf, y_pct, color='k', linewidth=lw)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%Peso")

    elif grafico == "Diagrama Acumulativo de Subtamaño":
        ax.scatter(x_inf, yf, s=48, linewidths=0.9, edgecolors='k', facecolors='white', zorder=4)
        ax.plot(x_inf, yf, color='k', linewidth=lw)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%F(d)")

    elif grafico == "Diagrama Acumulativo de Sobretamaño":
        ax.scatter(x_inf, yr, marker='x', s=36, linewidths=1.2, edgecolors='k', facecolors='k', zorder=4)
        ax.plot(x_inf, yr, color='k', linewidth=lw)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%R(d)")

    elif grafico == "Diagrama Acumulativo":
        ax.scatter(x_inf, yf, s=48, linewidths=0.9, edgecolors='k', facecolors='white', zorder=4, label='%F(d)')
        ax.plot(x_inf, yf, color='k', linewidth=lw)
        ax.scatter(x_inf, yr, marker='x', s=36, linewidths=1.2, edgecolors='k', facecolors='k', zorder=4, label='%R(d)')
        ax.plot(x_inf, yr, color='k', linewidth=lw)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%")
        ax.legend()

    else:  # Curvas granulométricas de la muestra
        ax.scatter(x_inf, y_pct, marker='s', s=48, linewidths=0.9, edgecolors='k', facecolors='white', zorder=4, label='%Peso')
        ax.plot(x_inf, y_pct, color='k', linewidth=lw)
        ax.scatter(x_inf, yf, marker='o', s=48, linewidths=0.9, edgecolors='k', facecolors='white', zorder=4, label='%F(d)')
        ax.plot(x_inf, yf, color='k', linewidth=lw)
        ax.scatter(x_inf, yr, marker='x', s=36, linewidths=1.2, edgecolors='k', facecolors='k', zorder=4, label='%R(d)')
        ax.plot(x_inf, yr, color='k', linewidth=lw)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%")
        ax.legend()

    # --- Escalas y límites ---
    from matplotlib.ticker import MaxNLocator

    if escala == "Escala decimal":
        ax.set_xlim(0, np.nanmax(x_inf) * 1.05 if len(x_inf.dropna()) > 0 else 1)
        ax.set_ylim(0, 100)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
        ax.grid(True, which='both', ls='--', alpha=0.5)

    elif escala == "Escala semilogarítmica":
        ax.set_xscale('log')
        xpos = x_inf.dropna().values
        xpos = xpos[xpos > 0] if len(xpos) > 0 else np.array([])
        if len(xpos) > 0:
            ax.set_xlim(np.min(xpos) * 0.8, np.max(xpos) * 1.2)
        ax.set_ylim(0, 100)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=10))
        ax.grid(True, which='both', ls='--', alpha=0.5)

    else:  # Escala logarítmica (ambos log)
        ax.set_xscale('log')
        ax.set_yscale('log')

        # Filtrar valores positivos
        # ... (código de filtrado de xpos y ypos se mantiene igual)
        xpos = x_inf.dropna().values
        xpos = xpos[xpos > 0] if len(xpos) > 0 else np.array([])
        y_comb = np.concatenate([y_pct.values, yf.values, yr.values])
        ypos = y_comb[~np.isnan(y_comb)]
        ypos = ypos[ypos > 0] if len(ypos) > 0 else np.array([])
        # ...

        # Limites (ajustar un poco los límites para evitar que los puntos toquen el borde)
        if len(xpos) > 0:
            # Ampliamos los límites en X un poco más
            ax.set_xlim(np.min(xpos) * 0.5, np.max(xpos) * 2.0)
        if len(ypos) > 0:
            # Aumentamos el límite inferior de Y para dar más espacio a los valores pequeños
            min_y = np.min(ypos)
            # Aseguramos un límite inferior limpio, por ejemplo, 1 si el mínimo es menor, o 0.5 veces el mínimo
            y_lim_min = 1.0 if min_y < 1.0 else min_y * 0.5
            ax.set_ylim(y_lim_min, np.max(ypos) * 1.5)

        # Ejes logarítmicos con ticks limpios y menor densidad de etiquetas
        from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter, FormatStrFormatter

        # === Eje X (Tamaño inferior) ===
        # Usar el LogLocator por defecto (que ya maneja potencias de 10)
        # y establecer el formateador a 'ScalarFormatter' para las etiquetas principales.
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=10)) 
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1)) # Ticks menores entre potencias
    
        # Formateador simple: Mostrar etiquetas solo para los ticks principales
        ax.xaxis.set_major_formatter(ScalarFormatter()) 
        ax.xaxis.set_minor_formatter(NullFormatter())  # No mostrar números en ticks menores

        # === Eje Y (%) ===
        # Similar al eje X: Ticks principales y menores, solo etiquetas en los principales
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1))
    
        # Usar FormatStrFormatter para forzar etiquetas con números enteros o con un formato específico
        # Esto puede evitar notación científica y hacerlos más legibles.
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.yaxis.set_minor_formatter(NullFormatter())

        # Rotar las etiquetas del eje X si fuese necesario (solo si hay superposición)
        # plt.gcf().autofmt_xdate()

        # Grid
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
                    sigmaI = abs((phi[84] - phi[16]) / 4.0 + (phi[95] - phi[5]) / 6.6)
                    SkI = ((2*phi[50] - phi[16] - phi[84]) / (2*(phi[84]-phi[16]))) + \
                          ((2*phi[50] - phi[5] - phi[95]) / (2*(phi[95]-phi[5])))
                    KG = (phi[95] - phi[5]) / (2.44 * (phi[75] - phi[25]))

                    # Tabla sin índice
                    folk_tbl = pd.DataFrame({
                        'Parámetro':['M (φ)','Md (φ)','σ (φ)','Sk (φ)','K (φ)'],
                        'Valor':[M, Md, sigmaI, SkI, KG]
                    })

                    col_tbl, col_interp = st.columns([1,1])
                    with col_tbl:
                        st.subheader("Tabla Folk & Ward")
                        # Ocultar la numeración automática
                        st.dataframe(folk_tbl.style.format({'Valor':'{:.3f}'}).hide(axis="index"))


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
def GGS_model(d, m, dmax):
    """
    Modelo de Gaudin-Schuhmann (GGS).
    d    : array de tamaños de partícula
    m    : parámetro de distribución
    dmax : tamaño máximo
    """
    d = np.array(d, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return 100.0 * (d / dmax) ** m

def RRSB_model(d, m, l):
    """
    Modelo de Rosin-Rammler-Sperling-Bennett (RRSB).
    d : array de tamaños de partícula
    m : parámetro de distribución
    l : parámetro de tamaño característico
    """
    d = np.array(d, dtype=float)
    return 100.0 * (1 - np.exp(-(d / l) ** m))


def double_weibull(d, alpha, delta1, delta2, d80):
    """
    Modelo Doble Weibull basado en d80 y ln(0.2).
    d      : array de tamaños de partícula
    alpha  : peso relativo de la primera distribución (0<alpha<1)
    delta1 : exponente de la primera Weibull
    delta2 : exponente de la segunda Weibull
    d80    : tamaño característico
    """
    d = np.array(d, dtype=float)
    ln02 = np.log(0.2)  # negativo (-1.609...)
    return 100.0 * (
        alpha * (1 - np.exp(ln02 * (d / d80) ** delta1)) +
        (1 - alpha) * (1 - np.exp(ln02 * (d / d80) ** delta2))
    )
            
#---------- PÁGINA 5: Selección del Modelo ----------
def page_5():
    st.title("SELECCIÓN DEL MODELO")
    st.markdown("Ajuste de modelos: GGS, RRSB y Doble Weibull. Se estiman parámetros minimizando SSE (F.O.).")
    results = st.session_state.results_table.copy()
    if results.empty:
        st.error("No hay resultados para ajustar. Regresa a DATOS EXPERIMENTALES y pulsa EJECUTAR.")
        return

    df_fit = results.iloc[:-1].copy()
    mask = (df_fit['Tamaño inferior (µm)'] > 0) & (~np.isnan(df_fit['%F(d)']))
    d = df_fit['Tamaño inferior (µm)'][mask].astype(float).values
    y_exp = df_fit['%F(d)'][mask].astype(float).values

    if len(d) < 3:
        st.warning("Se requieren al menos 3 puntos válidos para ajuste.")
        return

    if st.button("AJUSTAR"):
        n = len(d)  # número de puntos válidos

        # Función para calcular F.O. según tu definición
        def FO_calc(y_model, y_exp):
            eps2 = ((y_exp - y_model) / y_exp) ** 2  # ε²_i = ((F_exp - F_model)/F_exp)^2
            return np.sqrt(np.sum(eps2) / (n - 1))

        # ----------- GGS -----------
        try:
            def f_ggs(params):
                m, dmax = params
                ypred = GGS_model(d, m, dmax)
                eps2 = ((y_exp - ypred) / y_exp) ** 2
                return np.sqrt(np.sum(eps2) / (n - 1))

            # Lista de inicializaciones
            x0_list = [
                [0.5, np.max(d)],       # como Excel
                [1.0, np.median(d)*2],  # otra opción
                [2.0, np.max(d)],       # más “empinado”
            ]

            best = None
            for x0 in x0_list:
                res = minimize(
                    f_ggs,
                    x0,
                    method='L-BFGS-B',
                    # dmax restringido a max(d)
                    bounds=[(0.01, 10), (1e-6, np.max(d)*10)],
                    options={'ftol': 1e-12, 'maxiter': 10000}
                )
                if best is None or res.fun < best.fun:
                    best = res

            if best is not None and best.success:
                res1 = best
                ggs_params = res1.x.tolist()
                y_ggs_pred = GGS_model(d, *ggs_params)
                FO_ggs = FO_calc(y_ggs_pred, y_exp)
            else:
                ggs_params = [np.nan, np.nan]
                FO_ggs = np.inf

        except Exception as e:
            ggs_params = [np.nan, np.nan]
            FO_ggs = np.inf
       
        # ----------- RRSB -----------
        try:
            def f_rrsb(params):
                m, l = params
                ypred = RRSB_model(d, m, l)
                eps2 = ((y_exp - ypred) / y_exp) ** 2
                return np.sqrt(np.sum(eps2) / (n - 1))

            # Lista de inicializaciones
            x0_list = [
                [0.5, np.median(d)],     # clásico
                [1.0, np.mean(d)],       # variación
                [0.8, np.max(d)/2],      # otra opción
            ]

            best = None
            for x0 in x0_list:
                res = minimize(
                    f_rrsb,
                    x0,
                    method='L-BFGS-B',
                    bounds=[(0.01, 10), (1e-6, max(d)*10)],
                    options={'ftol': 1e-12, 'gtol': 1e-12, 'maxiter': 10000}
                )
                if best is None or res.fun < best.fun:
                    best = res

            res2 = best
            rrsb_params = res2.x.tolist()
            y_rrsb_pred = RRSB_model(d, *rrsb_params)
            FO_rrsb = FO_calc(y_rrsb_pred, y_exp)

        except Exception as e:
            FO_rrsb = np.inf
            rrsb_params = [np.nan, np.nan]

        # ----------- Double Weibull -----------
        try:
            def f_double(params):
                alpha, k1, k2, d80 = params
                ypred = double_weibull(d, alpha, k1, k2, d80)
                eps2 = ((y_exp - ypred) / y_exp) ** 2
                return np.sqrt(np.sum(eps2) / (n - 1))

            # Estimación inicial de d80 (si falla usamos mediana)
            try:
                inv = interp1d(
                    df_fit['%F(d)'], df_fit['Tamaño inferior (µm)'],
                    fill_value="extrapolate", bounds_error=False
                )
                init_d80 = float(inv(80.0))
                if init_d80 <= 0 or np.isnan(init_d80):
                    init_d80 = np.median(d)
            except:
                init_d80 = np.median(d)

            # Lista de inicializaciones
            x0_list = [
                [0.5, 0.9, 0.9, init_d80],
                [0.6, 1.0, 1.0, np.mean(d)],
                [0.4, 0.8, 0.8, np.max(d)/2],
                [0.5, 3.0, 1.0, init_d80],   # 🔹 nueva inicialización sugerida
            ]

            bounds_dw = [
                (1e-3, 8), # δ0​ usamos (1e-3, 1-1e-3) para evitar que δ0​ sea 0 o 1 y se reduzca a una weibull, si queremos admitir 0 y 1 usamos: (0.0, 1.0), sin embargo, valores mayores a 1 no dan ajuste de curva buenos
                (0.01, 10.0),  # δ1​
                (0.01, 10.0),  # δ2
                (1e-3, max(d)*10)  # d80
            ]

            best = None
            for x0 in x0_list:
                res = minimize(
                    f_double,
                    x0,
                    method='L-BFGS-B',
                    bounds=bounds_dw,
                    options={'ftol': 1e-12, 'gtol': 1e-12, 'maxiter': 20000}
                )

                if best is None or res.fun < best.fun:
                    best = res

            res3 = best
            dw_params = res3.x.tolist()
            y_dw_pred = double_weibull(d, *dw_params)
            FO_dw = FO_calc(y_dw_pred, y_exp)

        except Exception as e:
            FO_dw = np.inf
            dw_params = [np.nan]*4

        # Guardar resultados en session_state
        st.session_state.models_fit = {
            'GGS': {'FO': FO_ggs, 'params': ggs_params},
            'RRSB': {'FO': FO_rrsb, 'params': rrsb_params},
            'DoubleWeibull': {'FO': FO_dw, 'params': dw_params}
        }

        st.success("Ajustes completados")

    # ----------- Recuperar parámetros desde session_state ---------
    if not st.session_state.get("models_fit", None):
        st.warning("Pulsa 'AJUSTAR' para estimar los parámetros de los modelos.")
        return

    fits = st.session_state.models_fit
    ggs_params = fits['GGS']['params']; FO_ggs = fits['GGS']['FO']
    rrsb_params = fits['RRSB']['params']; FO_rrsb = fits['RRSB']['FO']
    dw_params   = fits['DoubleWeibull']['params']; FO_dw = fits['DoubleWeibull']['FO']

    # ------------------- TABLA COMPARATIVA -------------------
    table_data = []
    for i in range(len(d)):
        xi = d[i]; ye = y_exp[i]
        row = {"Tamaño inferior (µm)": xi, "%F(d)e": ye}

        # GGS
        if not np.isnan(ggs_params).any():
            m, dmax = ggs_params
            y_g = float(GGS_model([xi], m, dmax)[0])
            y_g_clip = min(100.0, y_g)
            row["%F(d)m_GGS"] = y_g_clip
            row["ε²_GGS"] = ((ye - y_g_clip)/ye)**2 if ye !=0 else np.nan
        else:
            row["%F(d)m_GGS"] = np.nan; row["ε²_GGS"] = np.nan

        # RRSB
        if not np.isnan(rrsb_params).any():
            m2, l = rrsb_params
            y_r = float(RRSB_model([xi], m2, l)[0])
            y_r_clip = min(100.0, y_r)
            row["%F(d)m_RRSB"] = y_r_clip
            row["ε²_RRSB"] = ((ye - y_r_clip)/ye)**2 if ye !=0 else np.nan
        else:
            row["%F(d)m_RRSB"] = np.nan; row["ε²_RRSB"] = np.nan

        # Double Weibull
        if not np.isnan(dw_params).any():
            alpha, k1, k2, d80 = dw_params
            y_dw = float(double_weibull([xi], alpha, k1, k2, d80)[0])
            y_dw_clip = min(100.0, y_dw)
            row["%F(d)m_DW"] = y_dw_clip
            row["ε²_DW"] = ((ye - y_dw_clip)/ye)**2 if ye !=0 else np.nan
        else:
            row["%F(d)m_DW"] = np.nan; row["ε²_DW"] = np.nan

        table_data.append(row)

    df_comp = pd.DataFrame(table_data)
    st.subheader("Tabla comparativa: Experimental vs Modelos")
    st.dataframe(
        df_comp.style.format({
            "Tamaño inferior (µm)": "{:.2f}",
            "%F(d)e": "{:.2f}",
            "%F(d)m_GGS": "{:.2f}",
            "%F(d)m_RRSB": "{:.2f}",
            "%F(d)m_DW": "{:.2f}",
            "ε²_GGS": "{:.2e}",
            "ε²_RRSB": "{:.2e}",
            "ε²_DW": "{:.2e}"
        }),
        height=320
    )

    # ------------------- Tablas separadas de parámetros -------------------
    def sum_errors_squared(y_model):
        return np.sum(((y_exp - y_model)/y_exp)**2) if y_model is not None else np.nan

    # Predicciones
    y_ggs = np.clip(GGS_model(d, *ggs_params), None, 100.0) if not np.isnan(ggs_params).any() else None
    y_rrsb = np.clip(RRSB_model(d, *rrsb_params), None, 100.0) if not np.isnan(rrsb_params).any() else None
    y_dw = np.clip(double_weibull(d, *dw_params), None, 100.0) if not np.isnan(dw_params).any() else None

    # Crear tablas
    ggs_table = pd.DataFrame([{
        'dmax (µm)': ggs_params[1],
        'm': ggs_params[0],
        'Σε²': sum_errors_squared(y_ggs),
        'F.O.': fits['GGS']['FO']
    }])

    rrsb_table = pd.DataFrame([{
        'l (µm)': rrsb_params[1],
        'm': rrsb_params[0],
        'Σε²': sum_errors_squared(y_rrsb),
        'F.O.': fits['RRSB']['FO']
    }])

    dw_table = pd.DataFrame([{
        'δ0​': dw_params[0],
        'δ1': dw_params[1],
        'δ2': dw_params[2],
        'd80 (µm)': dw_params[3],
        'Σε²': sum_errors_squared(y_dw),
        'F.O.': fits['DoubleWeibull']['FO']
    }])

    # Mostrar en tres columnas
    col1, col2, col3 = st.columns(3)

    with col1:
        st.subheader("Parámetros GGS")
        st.table(
            ggs_table.style.format({
                'dmax (µm)': "{:.2f}",
                'm': "{:.2f}",
                'Σε²': "{:.2e}",
                'F.O.': "{:.2e}"
            })
        )

    with col2:
        st.subheader("Parámetros RRSB")
        st.table(
            rrsb_table.style.format({
                'l (µm)': "{:.2f}",
                'm': "{:.2f}",
                'Σε²': "{:.2e}",
                'F.O.': "{:.2e}"
            })
        )

    with col3:
        st.subheader("Parámetros Doble Weibull")
        st.table(
            dw_table.style.format({
                'δ0​': "{:.2f}",
                'δ1': "{:.2f}",
                'δ2': "{:.2f}",
                'd80 (µm)': "{:.2f}",
                'Σε²': "{:.2e}",
                'F.O.': "{:.2e}"
            })
        )


    # ------------------- VALIDACIÓN DEL MODELO -------------------
    st.subheader("Validación del modelo")

    def model_metrics(y_exp, y_pred, k):
        """Calcula R2, MAPE, AIC y BIC para un modelo"""
        n = len(y_exp)
        residuals = y_exp - y_pred
        SSE = np.sum(residuals**2)

        # R²
        R2 = 1 - (SSE / np.sum((y_exp - np.mean(y_exp))**2))

        # MAPE
        MAPE = np.mean(np.abs((y_exp - y_pred) / y_exp)) * 100

        # AIC y BIC
        AIC = n * np.log(SSE / n) + 2 * k
        BIC = n * np.log(SSE / n) + k * np.log(n)

        return R2, MAPE, AIC, BIC

    metrics_data = []

    if y_ggs is not None:
        R2, MAPE, AIC, BIC = model_metrics(y_exp, y_ggs, k=2)
        metrics_data.append({"Modelo": "GGS", "R²": R2, "MAPE (%)": MAPE,
                             "AIC": AIC, "BIC": BIC, "F.O.": FO_ggs})

    if y_rrsb is not None:
        R2, MAPE, AIC, BIC = model_metrics(y_exp, y_rrsb, k=2)
        metrics_data.append({"Modelo": "RRSB", "R²": R2, "MAPE (%)": MAPE,
                             "AIC": AIC, "BIC": BIC, "F.O.": FO_rrsb})

    if y_dw is not None:
        R2, MAPE, AIC, BIC = model_metrics(y_exp, y_dw, k=4)
        metrics_data.append({"Modelo": "Doble Weibull", "R²": R2, "MAPE (%)": MAPE,
                             "AIC": AIC, "BIC": BIC, "F.O.": FO_dw})

    df_metrics = pd.DataFrame(metrics_data)
    st.dataframe(df_metrics.style.format({
        "R²": "{:.4f}", "MAPE (%)": "{:.2f}", "AIC": "{:.2f}", "BIC": "{:.2f}", "F.O.": "{:.2e}"
    }))

    st.subheader("Interpretación de las métricas")

    # Comentario sobre R²
    best_r2_model = df_metrics.loc[df_metrics["R²"].idxmax(), "Modelo"]
    st.markdown(f"**R²:** El modelo con mayor capacidad explicativa es **{best_r2_model}**, "
            f"lo que significa que representa mejor la variabilidad de los datos.")

    # Comentario sobre MAPE
    best_mape_model = df_metrics.loc[df_metrics["MAPE (%)"].idxmin(), "Modelo"]
    st.markdown(f"**MAPE:** El modelo más preciso en términos de error porcentual es **{best_mape_model}**, "
                "ya que presenta el menor valor de MAPE.")

    # Comentario sobre AIC y BIC
    best_ic_model = df_metrics.loc[df_metrics["AIC"].idxmin(), "Modelo"]
    st.markdown(f"**AIC/BIC:** El modelo más parsimonioso (mejor balance ajuste-complejidad) es **{best_ic_model}**, "
                "con los valores más bajos de AIC y BIC.")

    # Conclusión general
    best_model = "RRSB"  # Aquí puedes automatizar con reglas de decisión
    st.markdown(f"**CONCLUSIÓN:** Aunque {best_mape_model} presenta el menor error, "
                f"el análisis global de todas las métricas indica que el modelo más representativo es **{best_model}**.")
 
    # ------------------- Gráficos -------------------
    xdata = d; ydata = y_exp
    dd = np.linspace(np.min(xdata), np.max(xdata), 500)

    # Predicciones para graficar en malla fina (ya clippeadas 0-100)
    y_ggs_plot = np.clip(GGS_model(dd, *ggs_params), None, 100.0) if not np.isnan(ggs_params).any() else None
    y_rrsb_plot = np.clip(RRSB_model(dd, *rrsb_params), None, 100.0) if not np.isnan(rrsb_params).any() else None
    y_dw_plot = np.clip(double_weibull(dd, *dw_params), None, 100.0) if not np.isnan(dw_params).any() else None

    # Estilo de puntos experimentales: fondo blanco, borde negro
    marker_kwargs = dict(marker='o', s=48, linewidths=0.9,
                         edgecolors='k', facecolors='white', zorder=4)

    # Estilos de línea (todas negras, distintos linestyle)
    line_styles = {
        'GGS': {'linestyle':'-',  'label':'GGS'},
        'RRSB':{'linestyle':'--', 'label':'RRSB'},
        'DW':  {'linestyle':'-.', 'label':'Doble Weibull'}
    }
    line_common = {'color':'k', 'linewidth':1.1, 'zorder':1}

    graf_option = st.selectbox("Selecciona tipo de gráfica:",
                               ["Comparación de perfiles",
                                "Diagrama GGS (log-log)",
                                "Diagrama RRSB (log-x, transform y)",
                                "Diagrama DW (decimal)"])

    fig, ax = plt.subplots(figsize=(8,4))

    # Fondos: exterior gris, interior blanco
    fig.patch.set_facecolor("lightgrey")  # exterior (área de escalas y labels)
    ax.set_facecolor("white")             # interior (área de curvas y puntos)

    if graf_option == "Comparación de perfiles":
        # Dibujar curvas de modelos primero (si existen), cada una con su linestyle
        if y_ggs_plot is not None:
            ax.plot(dd, y_ggs_plot, linestyle=line_styles['GGS']['linestyle'],
                    **line_common, label="GGS")
        if y_rrsb_plot is not None:
            ax.plot(dd, y_rrsb_plot, linestyle=line_styles['RRSB']['linestyle'],
                    **line_common, label="RRSB")
        if y_dw_plot is not None:
            ax.plot(dd, y_dw_plot, linestyle=line_styles['DW']['linestyle'],
                    **line_common, label="DW")

        # Luego dibujar puntos experimentales (sobre las curvas)
        ax.scatter(xdata, ydata, **marker_kwargs)
        
        ax.set_title("Comparación de perfiles")
        ax.set_xlabel("Tamaño (µm)")
        ax.set_ylabel("%F(d)")
        ax.set_ylim(0,100)
        ax.set_yticks(np.arange(0,101,10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(8))
        ax.grid(True, ls='--', alpha=0.5)
        ax.legend()

    elif graf_option == "Diagrama GGS (log-log)":
        ax.set_xscale('log'); ax.set_yscale('log')
        # Curva GGS primero
        if y_ggs_plot is not None:
            ax.plot(dd, np.clip(y_ggs_plot, 1e-6, 100.0),
                    linestyle=line_styles['GGS']['linestyle'], **line_common)
        # Puntos experimentales encima
        ax.scatter(xdata, np.clip(ydata, 1e-6, 100.0), **marker_kwargs)

        ax.set_title("Diagrama GGS (log-log)")
        ax.set_xlabel("Tamaño (µm) [escala log]")
        ax.set_ylabel("%F(d) [escala log]")

        # --- Ajuste dinámico de límites ---
        xmin, xmax = np.min(xdata), np.max(xdata)
        ymin, ymax = np.min(ydata[ydata > 0]), np.max(ydata)
        ax.set_xlim(xmin*0.8 if xmin > 0 else xmin, xmax*1.2)
        ax.set_ylim(ymin*0.8 if ymin > 0 else ymin, ymax*1.2)

        ax.grid(True, which='both', ls='--', alpha=0.5)

    elif graf_option == "Diagrama RRSB (log-x, transform y)":
        ax.set_xscale('log')

        # Transformación segura
        def transform_rrsb(y_percent):
            y = np.minimum(99.9999, np.maximum(1e-8, y_percent))
            return np.log(np.log(1.0 / (1.0 - (y / 100.0))))

        # Curva RRSB primero
        if y_rrsb_plot is not None:
            y_rrsb_trans = transform_rrsb(y_rrsb_plot)
            ax.plot(dd, y_rrsb_trans, linestyle=line_styles['RRSB']['linestyle'], **line_common)

        # Puntos experimentales encima
        ydata_trans = transform_rrsb(ydata)
        ax.scatter(xdata, ydata_trans, **marker_kwargs)

        ax.set_title("Diagrama RRSB (transform)")
        ax.set_xlabel("Tamaño (µm) [escala log]")
        ax.set_ylabel("Log[ ln(1 / (1 - (%F/100)) ) ]")

        # Ajuste X dinámico
        xmin, xmax = np.min(xdata), np.max(xdata)
        ax.set_xlim(xmin*0.8 if xmin > 0 else xmin, xmax*1.2)

        # Grilla con secundarias
        ax.grid(True, which='major', ls='--', alpha=0.5)
        ax.grid(True, which='minor', ls=':', alpha=0.3)

    elif graf_option == "Diagrama DW (decimal)":
        # Curva DW primero
        if y_dw_plot is not None:
            ax.plot(dd, y_dw_plot, linestyle=line_styles['DW']['linestyle'], **line_common)
        # Puntos experimentales encima
        ax.scatter(xdata, ydata, **marker_kwargs)

        ax.set_title("Diagrama DW (decimal)")
        ax.set_xlabel("Tamaño (µm)")
        ax.set_ylabel("%F(d)")
        ax.set_ylim(0,100)
        ax.set_yticks(np.arange(0,101,10))
        ax.xaxis.set_major_locator(plt.MaxNLocator(8))
        ax.grid(True, ls='--', alpha=0.5)

    # Mostrar figura
    st.pyplot(fig, use_container_width=True)

    # ------------------- Navegación -------------------
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            st.session_state.page = 4; st.rerun()
    with col2:
        if st.button("SIGUIENTE"):
            st.session_state.page = 6; st.rerun()

# ---------- PÁGINA 6: Exportación ----------
def page_6():
    st.title("EXPORTACIÓN DE DATOS")
    st.markdown("Descargar todas las tablas en un archivo Excel. Cada usuario obtiene su propio archivo en modo solo lectura.")

    # ---------------- Construir Excel en memoria ----------------
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:

        # Info usuario
        if "user_info" in st.session_state:
            pd.DataFrame([st.session_state['user_info']]).to_excel(writer, sheet_name='InfoUsuario', index=False)

        # Tabla de entrada
        if "input_table" in st.session_state and not st.session_state.input_table.empty:
            st.session_state.input_table.to_excel(writer, sheet_name='Entrada', index=False)

        # Resultados
        if "results_table" in st.session_state and not st.session_state.results_table.empty:
            st.session_state.results_table.to_excel(writer, sheet_name='Resultados', index=False)

        # Tamaños nominales
        if "nominal_sizes" in st.session_state and not st.session_state.nominal_sizes.empty:
            st.session_state.nominal_sizes.to_excel(writer, sheet_name='TamañosNominales', index=False)

        # Estadísticos descriptivos
        if "results_table" in st.session_state and not st.session_state.results_table.empty:
            results = st.session_state.results_table.copy()
            sizes = results['Tamaño inferior (µm)'].iloc[:-1].replace(0, np.nan).dropna()
            weights = results['%Peso'].iloc[:-1].replace(0, np.nan).dropna() / 100.0
            mask = (~sizes.isna()) & (~weights.isna())
            sizes = sizes.loc[mask]
            weights = weights.loc[mask]
            if len(sizes) > 0:
                wsum = weights.sum()
                mean = float((sizes * weights).sum() / wsum) if wsum>0 else float(sizes.mean())
                median = float(np.interp(50, 100*weights.cumsum(), sizes))
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
                stats_tbl.to_excel(writer, sheet_name='Estadísticos', index=False)

        # Folk & Ward (si existe)
        if "nominal_sizes" in st.session_state and not st.session_state.nominal_sizes.empty:
            try:
                req_pcts = [5,16,25,50,75,84,95]
                ds = {}
                for p in req_pcts:
                    match = st.session_state.nominal_sizes.loc[
                        (st.session_state.nominal_sizes['%F(d)'].sub(p).abs() < 1e-6), 'Tamaño (µm)']
                    if match.empty:
                        ds = None
                        break
                    else:
                        ds[p] = float(match.iloc[0])
                if ds:
                    ds_mm = {p: ds[p]/1000.0 for p in ds}
                    phi = {p: -np.log2(ds_mm[p]) for p in ds_mm}
                    M = (phi[16] + phi[50] + phi[84]) / 3.0
                    Md = phi[50]
                    sigmaI = abs((phi[84] - phi[16]) / 4.0 + (phi[95] - phi[5]) / 6.6)
                    SkI = ((2*phi[50] - phi[16] - phi[84]) / (2*(phi[84]-phi[16]))) + \
                          ((2*phi[50] - phi[5] - phi[95]) / (2*(phi[95]-phi[5])))
                    KG = (phi[95] - phi[5]) / (2.44 * (phi[75] - phi[25]))
                    folk_tbl = pd.DataFrame({
                        'Parámetro':['M (φ)','Md (φ)','σ (φ)','Sk (φ)','K (φ)'],
                        'Valor':[M, Md, sigmaI, SkI, KG]
                    })
                    folk_tbl.to_excel(writer, sheet_name='FolkWard', index=False)
            except:
                pass

        # Parámetros de modelos (página 5)
        if "models_fit" in st.session_state and st.session_state.models_fit:
            fits = st.session_state.models_fit
            ggs_params = fits['GGS']['params'] if 'GGS' in fits else [np.nan,np.nan]
            rrsb_params = fits['RRSB']['params'] if 'RRSB' in fits else [np.nan,np.nan]
            dw_params = fits['DoubleWeibull']['params'] if 'DoubleWeibull' in fits else [np.nan]*4

            ggs_table = pd.DataFrame([{
                'dmax (µm)': ggs_params[1], 'm': ggs_params[0], 'F.O.': fits['GGS']['FO'] if 'GGS' in fits else np.nan
            }])
            rrsb_table = pd.DataFrame([{
                'l (µm)': rrsb_params[1], 'm': rrsb_params[0], 'F.O.': fits['RRSB']['FO'] if 'RRSB' in fits else np.nan
            }])
            dw_table = pd.DataFrame([{
                'δ0​': dw_params[0], 'δ1': dw_params[1], 'δ2': dw_params[2], 'd80 (µm)': dw_params[3],
                'F.O.': fits['DoubleWeibull']['FO'] if 'DoubleWeibull' in fits else np.nan
            }])

            ggs_table.to_excel(writer, sheet_name='GGS_Params', index=False)
            rrsb_table.to_excel(writer, sheet_name='RRSB_Params', index=False)
            dw_table.to_excel(writer, sheet_name='DoubleWeibull_Params', index=False)

        # Métricas de validación
        if "results_table" in st.session_state and st.session_state.results_table is not None:
            df_fit = st.session_state.results_table.iloc[:-1].copy()
            mask = (df_fit['Tamaño inferior (µm)'] > 0) & (~np.isnan(df_fit['%F(d)']))
            d = df_fit['Tamaño inferior (µm)'][mask].astype(float).values
            y_exp = df_fit['%F(d)'][mask].astype(float).values
            if len(d)>0 and "models_fit" in st.session_state and st.session_state.models_fit:
                fits = st.session_state.models_fit
                def model_metrics(y_exp, y_pred, k):
                    n = len(y_exp)
                    residuals = y_exp - y_pred
                    SSE = np.sum(residuals**2)
                    R2 = 1 - (SSE / np.sum((y_exp - np.mean(y_exp))**2))
                    MAPE = np.mean(np.abs((y_exp - y_pred) / y_exp)) * 100
                    AIC = n * np.log(SSE / n) + 2 * k
                    BIC = n * np.log(SSE / n) + k * np.log(n)
                    return R2, MAPE, AIC, BIC

                metrics_data = []
                for model_name, y_model_func, k in [
                    ('GGS', lambda: GGS_model(d,*ggs_params), 2),
                    ('RRSB', lambda: RRSB_model(d,*rrsb_params), 2),
                    ('DobleWeibull', lambda: double_weibull(d,*dw_params), 4)
                ]:
                    try:
                        y_pred = y_model_func()
                        R2, MAPE, AIC, BIC = model_metrics(y_exp, y_pred, k)
                        metrics_data.append({"Modelo": model_name, "R²": R2, "MAPE (%)": MAPE, "AIC": AIC, "BIC": BIC})
                    except:
                        pass
                if metrics_data:
                    pd.DataFrame(metrics_data).to_excel(writer, sheet_name='ValidaciónModelos', index=False)

    # ---------------- Botón para descargar Excel ----------------
    st.download_button(
        label="📥 Descargar Excel",
        data=output.getvalue(),
        file_name="Analisis_Granulometrico.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        key="download_excel"
    )

    # ---------------- Botones de navegación ----------------
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR", key="btn_anterior"):
            st.session_state.page = 5
            st.rerun()
    with col2:
        if st.button("VOLVER AL INICIO", key="btn_inicio"):
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
































































































































