# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 21:15:05 2025

@author: ALEX
"""

# app_streamlit_granulometria.py
# Streamlit app: CARACTERIZACIÓN GRANULOMÉTRICA
# Autor: Alex Quispe (estructura solicitada por el usuario)
# Instrucciones: subir este archivo a GitHub y publicar en Streamlit Community Cloud o ejecutar localmente:
# streamlit run app_streamlit_granulometria.py

import streamlit as st
import pandas as pd
import numpy as np
import io
import base64
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import qrcode
from PIL import Image
import tempfile
import openpyxl

st.set_page_config(page_title="CARACTERIZACIÓN GRANULOMÉTRICA", layout="wide")

# ---------- Datos fijos: Serie Tyler (malla -> apertura µm) ----------
TYLER_SERIES = {
    "1.05\"": None,  # placeholder header-like
}

# More complete mapping extracted and corrected from your prompt (keys: malla número, value: µm)
TYLER = {
    26:500,22:400,19:300,16:16000,13:13200,11:11200,9:9500,
    8:8000,7:6700,6:5600,5:4760,4:4000,3.5:3350,3:2800,2.5:2360,
    2:2000,1.7:1700,1.4:1400,1.18:1180,1.0:1000,0.85:850,0.71:710,
    0.6:600,0.5:500,0.425:425,0.355:355,0.3:300,0.25:250,0.212:212,
    0.18:180,0.15:150,0.125:125,0.106:106,0.09:90,0.075:75,0.063:63,
    0.053:53,0.045:45,0.038:38
}

# ---------- NUEVA FUNCIÓN DE ESTILO ----------
def aplicar_estilo_granulometrico(ax):
    """Aplica un estilo consistente a los gráficos granulométricos."""
    # Colores
    LINE_COLOR = '#004488'
    GRID_COLOR = 'gray'
    BG_COLOR = '#F0F2F6'  # Streamlit background
    
    # Fondo
    ax.set_facecolor(BG_COLOR)
    
    # Líneas de la cuadrícula
    ax.grid(True, which='both', linestyle='-', linewidth=0.5, color=GRID_COLOR, alpha=0.6)
    
    # Líneas y marcadores por defecto
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 8
    
    # Color de las espinas (bordes del gráfico)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(GRID_COLOR)
    ax.spines['bottom'].set_color(GRID_COLOR)

# ---------- Utilidades ----------
def size_label_from_malla(malla_number):
    """Return readable label for a malla key (which may be non-integer)."""
    return f"{malla_number}# ({TYLER.get(malla_number,'?')} µm)"

def ensure_session():
    if 'page' not in st.session_state:
        st.session_state.page = 1
    if 'user_info' not in st.session_state:
        st.session_state.user_info = {}
    if 'input_table' not in st.session_state:
        st.session_state.input_table = pd.DataFrame()
    if 'selected_mode' not in st.session_state:
        st.session_state.selected_mode = None
    if 'results_table' not in st.session_state:
        st.session_state.results_table = pd.DataFrame()
    if 'nominal_sizes' not in st.session_state:
        st.session_state.nominal_sizes = pd.DataFrame(columns=['%F(d)','Tamaño (µm)'])
    if 'models_fit' not in st.session_state:
        st.session_state.models_fit = {}
    if 'uploaded_image' not in st.session_state:
        st.session_state.uploaded_image = None

ensure_session()

# ---------- PÁGINA 1: Bienvenida ----------
def page_1():
    st.title("CARACTERIZACIÓN GRANULOMÉTRICA")
    st.markdown(f"**Desarrollado por:** Alex Quispe — UTO, Carrera de Metalurgia y Ciencia de Materiales.")
    st.markdown("""
    Este programa permite realizar un análisis granulométrico completo a partir de datos experimentales
    (pesos retenidos por tamiz o mediciones manuales).
    Genera tablas, gráficos en distintas escalas,
    estadísticas descriptivas e inferenciales básicas, ajuste a modelos (GGS, RRSB, Doble Weibull),
    y exportación de resultados.
    """)
    st.markdown("**Imagen representativa** (subir desde tu PC):")
    img_file = st.file_uploader("Selecciona una imagen (jpg, png). Se mostrará en la portada.", type=['png','jpg','jpeg'])
    if img_file is not None:
        img = Image.open(img_file)
        st.image(img, width=300)
        st.session_state.uploaded_image = img_file
    st.subheader("Información del usuario y muestra")
    with st.form("user_info_form"):
        nombre = st.text_input("USUARIO (tu nombre)", value=st.session_state.user_info.get('nombre',''))
        correo = st.text_input("CORREO (Gmail)", 
value=st.session_state.user_info.get('correo',''))
        procedencia = st.text_input("PROCEDENCIA DE LA MUESTRA", value=st.session_state.user_info.get('procedencia',''))
        codigo = st.text_input("CÓDIGO (muestra)", value=st.session_state.user_info.get('codigo',''))
        fecha_muestreo = st.date_input("FECHA DE MUESTREO", value=st.session_state.user_info.get('fecha', datetime.today().date()))
        inicio = st.form_submit_button("INICIO")
        if inicio:
            st.session_state.user_info.update({
                'nombre': nombre, 'correo': correo, 'procedencia': procedencia,
                'codigo': codigo, 'fecha': fecha_muestreo.isoformat()
            })
            st.session_state.page = 2
            st.rerun()

# ---------- PÁGINA 2: Datos experimentales ----------
def page_2():
    st.title("DATOS EXPERIMENTALES")
    st.markdown("Inserte tamaños y pesos retenidos sobre cada tamiz para efectuar el análisis granulométrico.")
    total_weight = st.number_input("Peso total (g)", min_value=0.0, value=1000.0, step=0.1, key='peso_total')
    st.markdown("Elija método para introducir datos:")
    mode = st.radio("Modo", ["SELECCIONAR MALLAS", "INSERTAR MANUALMENTE"], index=0)
    st.session_state.selected_mode = mode

    if mode == "SELECCIONAR MALLAS":
        st.info("Seleccione las mallas de la serie Tyler. Se generará una tabla con la abertura (µm) prellenada.")
        # Build options from TYLER dict, order descending by aperture
        malla_options = sorted(TYLER.items(), key=lambda x: -x[1])
        labels = [f"{int(k) if float(k).is_integer() else k} - {v} µm" for k,v in malla_options]
        # Use multiselect
        selected = st.multiselect("Selecciona mallas (múltiple)", labels)
        # Map selected labels back to keys
        selected_keys = []
        for lab in selected:
            key = lab.split(" - ")[0]
            try:
                k = int(key)
            except:
                try:
                    k = float(key)
                except:
                    k = key
            selected_keys.append(k)
        # Generate default table
        if 'generated_mallas' not in st.session_state:
            st.session_state.generated_mallas = []
        if st.button("Generar tabla de mallas"):
            rows = []
            for k in selected_keys:
                rows.append({'Nº Malla (Tyler)': str(k)+'#', 'Abertura (µm)': TYLER.get(k, np.nan), 'Peso (g)': np.nan})
            df = pd.DataFrame(rows)
            st.session_state.input_table = df
            st.session_state.generated_mallas = selected_keys
            st.rerun()

    else:  # INSERTAR MANUALMENTE
        st.info("Ingrese el número de datos (entre 3 y 25).")
        n = st.number_input("Número de filas a insertar", min_value=3, max_value=25, value=6, step=1)
        if st.button("Generar tabla manual"):
            df = pd.DataFrame({
                'Tamaño (µm)': [np.nan]*n,
                'Peso (g)': [np.nan]*n
            })
            st.session_state.input_table = df
            st.rerun()

    st.markdown("**Tabla de entrada** (edítala con los pesos y tamaños necesarios):")
    # Show editable 
    if not st.session_state.input_table.empty:
        edited = st.data_editor(st.session_state.input_table, num_rows="dynamic")
        st.session_state.input_table = edited

    st.markdown("**Recomendación:** seleccionar MALLAS si trabajas con tamices; insertar manualmente si no se usan mallas estandarizadas.")

    if st.button("EJECUTAR"):
        # Validate
        df = st.session_state.input_table.copy()
        if df.empty:
            st.error("Debes generar y completar la tabla antes de ejecutar.")
        else:
            # Ensure Peso column exists
            if 'Peso (g)' not in df.columns and 'Peso (g)' in df.columns:
                pass
            st.session_state.page = 3
            st.rerun()

# ---------- Helper: compute granulometric analysis ----------
def compute_analysis(df_in, mode, total_weight):
    """
    df_in: DataFrame con columnas ['Nº Malla (Tyler)','Abertura (µm)','Peso (g)']
           o ['Tamaño (µm)','Peso (g)']
    total_weight: valor ingresado en page_2
    """
    df = df_in.copy()

    # Normalizar nombres
    if 'Abertura (µm)' in df.columns:
        df.rename(columns={'Abertura (µm)': 'Tamaño inferior (µm)',
                           'Peso (g)': 'Peso (g)'}, inplace=True)
    if 'Tamaño (µm)' in df.columns and 'Tamaño inferior (µm)' not in df.columns:
        df.rename(columns={'Tamaño (µm)': 'Tamaño inferior (µm)'}, inplace=True)

    # Asegurar datos numéricos
    df['Tamaño inferior (µm)'] = pd.to_numeric(df['Tamaño inferior (µm)'], errors='coerce')
    df['Peso (g)'] = pd.to_numeric(df['Peso (g)'], errors='coerce').fillna(0.0)

    # Ordenar de mayor a menor tamaño
    df = df.sort_values(by='Tamaño inferior (µm)', ascending=False).reset_index(drop=True)

    # Calcular Tamaño superior
    size_sup = []
    for i in range(len(df)):
        if i == 0:
            sup = df.loc[i, 'Tamaño inferior (µm)'] * np.sqrt(2)
        else:
            sup = df.loc[i-1, 'Tamaño inferior (µm)']
        size_sup.append(sup)
    df['Tamaño superior (µm)'] = size_sup

    # ---- Fila extra_row (< Tamaño mínimo) ----
    peso_resto = max(total_weight - df['Peso (g)'].sum(), 0.0)
    if len(df) > 0:
        extra_row = {
            'Tamaño inferior (µm)': 0.0,
            'Tamaño superior (µm)': df.loc[len(df)-1, 'Tamaño inferior (µm)'],
            'Tamaño promedio (µm)': df.loc[len(df)-1, 'Tamaño inferior (µm)'] / 2,
            'Peso (g)': peso_resto
        }
        df = pd.concat([df, pd.DataFrame([extra_row])], ignore_index=True)

    # Tamaño promedio
    df['Tamaño promedio (µm)'] = (df['Tamaño superior (µm)'] + df['Tamaño inferior (µm)']) / 2.0

    # %Peso respecto al total ingresado
    df['%Peso'] = 100.0 * df['Peso (g)'] / total_weight

    # %R(d) acumulado retenido y %F(d) pasante
    df['%R(d)'] = df['%Peso'].cumsum()
    df['%F(d)'] = 100.0 - df['%R(d)']

    # ---- Fila de TOTALES ----
    total_row = {
        'Tamaño superior (µm)': np.nan,
        'Tamaño inferior (µm)': np.nan,
        'Tamaño promedio (µm)': np.nan,
        'Peso (g)': total_weight,
        '%Peso': 100.0,
        '%R(d)': np.nan,
        '%F(d)': np.nan
    }
    if mode == "SELECCIONAR MALLAS" and 'Nº Malla (Tyler)' in df_in.columns:
        df['Nº de malla (intervalo)'] = list(df_in['Nº Malla (Tyler)']) + ['extra']  # para el resto
        total_row['Nº de malla (intervalo)'] = 'TOTAL'
        cols_order = ['Nº de malla (intervalo)', 'Tamaño superior (µm)',
                      'Tamaño inferior (µm)', 'Tamaño promedio (µm)',
                      'Peso (g)', '%Peso', '%F(d)', '%R(d)']
    else:
        cols_order = ['Tamaño superior (µm)', 'Tamaño inferior (µm)',
                      'Tamaño promedio (µm)', 'Peso (g)',
                      '%Peso', '%F(d)', '%R(d)']

    # Concatenar fila de totales
    df = pd.concat([df, pd.DataFrame([total_row])], ignore_index=True)

    return df[cols_order]

# ---------- PÁGINA 3: Análisis granulométrico (Resultados) ----------
def page_3():
    st.title("ANÁLISIS GRANULOMÉTRICO")
    df_in = st.session_state.input_table.copy()
    total_weight = st.session_state.get('peso_total', 0.0)

    if df_in.empty:
        st.error("No hay datos de entrada. Regresa y genera la tabla de datos.")
        if st.button("Regresar"):
            st.session_state.page = 2
            st.rerun()
        return

    # Evitar que se recalculen los resultados cada vez que cambias gráfico
    if st.session_state.results_table.empty:
        results = compute_analysis(df_in, st.session_state.selected_mode, total_weight)
        st.session_state.results_table = results
    else:
        results = st.session_state.results_table

    st.markdown("**Tabla de resultados**")
    st.dataframe(
        results.style.format({
            "Peso (g)": "{:.3f}",
            "%Peso": "{:.3f}",
            "%F(d)": "{:.3f}",
            "%R(d)": "{:.3f}"
        }),
        height=300
    )

    # Gráficos
    st.markdown("**Seleccione gráfico**")
    grafico = st.selectbox("SELECCIONE GRÁFICO", [
        "Histograma de frecuencia",
        "Diagrama de simple distribución",
        "Diagrama Acumulativo de Subtamaño",
        "Diagrama Acumulativo de Sobretamaño",
        "Diagrama Acumulativo (Combinación)",
        "Curvas granulométricas (Combinación 2,3,4)"
    ])
    escala = st.selectbox(
        "Escala",
        ["Escala decimal", "Escala semilogarítmica (X log)", "Escala logarítmica (ambos log)"]
    )

    # Datos para gráficos (AHORA SE USA 'Tamaño inferior (µm)')
    plot_df = results.copy()
    x = plot_df['Tamaño inferior (µm)'].replace(0, np.nan).dropna()
    y_pct = plot_df['%Peso'].loc[x.index]
    yf = plot_df['%F(d)'].loc[x.index]
    yr = plot_df['%R(d)'].loc[x.index]
    
    fig, ax = plt.subplots(figsize=(8, 4))
    aplicar_estilo_granulometrico(ax)

    # Plot based on selected graph type
    if grafico == "Histograma de frecuencia":
        ax.bar(x, y_pct, width=np.nanmax(x) / len(x) if len(x) > 0 else 1, color='#004488')
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%Peso")
    elif grafico == "Diagrama de simple distribución":
        ax.plot(x, y_pct, marker='s', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%Peso")
    elif grafico == "Diagrama Acumulativo de Subtamaño":
        ax.plot(x, yf, marker='o', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%F(d)")
    elif grafico == "Diagrama Acumulativo de Sobretamaño":
        ax.plot(x, yr, marker='x', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("%R(d)")
    elif grafico == "Diagrama Acumulativo (Combinación)":
        ax.plot(x, yf, label='%F(d)', marker='o', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.plot(x, yr, label='%R(d)', marker='x', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("Porcentaje")
        ax.legend()
    else:  # Curvas granulométricas (2,3,4)
        ax.plot(x, y_pct, label='%Peso', marker='s', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.plot(x, yf, label='%F(d)', marker='o', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.plot(x, yr, label='%R(d)', marker='x', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.set_xlabel("Tamaño inferior (µm)")
        ax.set_ylabel("Porcentaje")
        ax.legend()

    # Escalas
    if escala == "Escala semilogarítmica (X log)":
        ax.set_xscale('log')
    elif escala == "Escala logarítmica (ambos log)":
        ax.set_xscale('log')
        ax.set_yscale('log')
        
    ax.set_title(grafico)
    st.pyplot(fig)

    if st.button("SIGUIENTE"):
        st.session_state.page = 4
        st.rerun()

# ---------- PÁGINA 4: Análisis de datos ----------
def page_4():
    st.title("ANÁLISIS DE DATOS")
    results = st.session_state.results_table.copy()
    if results.empty:
        st.error("No hay resultados calculados.")
        if st.button("Regresar"):
            st.session_state.page = 3
            st.rerun()
        return

    # ---------- Sección 1: ANÁLISIS ESTADÍSTICO ----------
    st.header("ANÁLISIS ESTADÍSTICO")
    # Use Tamaño promedio weighted by %Peso to compute distribution statistics
    sizes = results['Tamaño promedio (µm)'].replace(0, np.nan).dropna()
    weights = results['%Peso'].loc[sizes.index] / 100.0  # fraction
    # If weights sum to zero, fallback to equal weight
    if weights.sum() <= 0:
        weights = np.ones_like(sizes) / len(sizes)
    mean = np.sum(sizes * weights) / np.sum(weights)
    median = np.nanpercentile(np.repeat(sizes.values, (weights*1000).astype(int)+1), 50) if len(sizes)>0 else np.nan
    # approximate mode: size with max %Peso
    mode = sizes.iloc[np.argmax(results['%Peso'].values[:len(sizes)])] if len(sizes)>0 else np.nan
    variance = np.sum(weights * (sizes - mean)**2) / np.sum(weights)
    std = np.sqrt(variance)
    kurtosis = pd.Series(sizes).kurtosis()
    rango = np.nanmax(sizes) - np.nanmin(sizes) if len(sizes)>0 else np.nan
    stats_tbl = pd.DataFrame({
        'Estadístico':['Media (µm)','Mediana (µm)','Moda (µm)','Varianza (µm²)','Desvío estándar (µm)','Curtosis','Rango (µm)'],
        'Valor':[mean, median, mode, 
variance, std, kurtosis, rango]
    })
    st.table(stats_tbl)

    st.markdown("**Interpretación**: La varianza indica la dispersión de tamaños; un valor mayor significa mayor dispersión. El rango (máximo - mínimo) muestra la amplitud de tamaños presentes en la muestra.")

    # ---------- Sección 2: REGISTRO DE TAMAÑOS Y PORCENTAJES ----------
    st.header("REGISTRO DE TAMAÑOS Y PORCENTAJES")
    cols = st.columns(2)
    with cols[0]:
        st.subheader("PERFIL GRANULOMÉTRICO")
        st.markdown("El perfil granulométrico muestra la curva acumulativa de subtamaño (%F(d)).")
        
        fig, ax = plt.subplots(figsize=(5,4))
        aplicar_estilo_granulometrico(ax)
        
        x = results['Tamaño promedio (µm)'].replace(0, np.nan).dropna()
        y = results['%F(d)'].loc[x.index]
        ax.plot(x, y, marker='o', ls='-', mfc='white', mec='#004488', color='#004488', lw=2, ms=8)
        ax.set_xlabel("Tamaño (µm)")
        ax.set_ylabel("%F(d)")
        ax.set_xscale('linear')
        ax.set_yscale('linear')
        ax.grid(True)
        st.pyplot(fig)
        
        # compute d80 by interpolation/extrapolation
        try:
            mask = ~np.isnan(x) & ~np.isnan(y)
            xs = np.array(x[mask])
            ys = np.array(y[mask])
            
            if len(xs) >= 2:
                order = np.argsort(xs)
                xs_s = xs[order]
                ys_s = ys[order]
                
                if (ys_s.min() <= 80 <= ys_s.max()) or (ys_s.max() >= 80 and ys_s.min() <= 80):
                    inv = interp1d(ys_s, xs_s, fill_value="extrapolate")
                    d80 = float(inv(80.0))
                    st.markdown(f"**d80 = {d80:.3f} µm**")
                    # Flecha eliminada
                else:
                    inv = interp1d(ys_s, xs_s, fill_value="extrapolate")
                    d80 = float(inv(80.0))
                    st.warning("Los datos no cubren 80% — se realizó EXTRAPOLACIÓN.")
                    st.markdown(f"**d80 (extrapolado) = {d80:.3f} µm**")
            else:
                st.warning("No hay suficientes datos para calcular d80.")
        except Exception as e:
            st.error("Error calculando d80: " + str(e))

    with cols[1]:
        st.subheader("TAMAÑOS NOMINALES")
        if st.button("CALCULAR TAMAÑOS NOMINALES"):
            pct = st.number_input("Ingrese %F(d) =", min_value=0.0, max_value=100.0, value=50.0, key='calc_nom_pct')
            if st.button("GRABAR %->Tamaño", key='grabar_nom'):
                x = results['Tamaño promedio (µm)'].replace(0, np.nan)
                y = results['%F(d)']
                mask = ~np.isnan(x) & ~np.isnan(y)
                try:
                    inv = interp1d(y[mask], x[mask], fill_value="extrapolate")
                    d_val = float(inv(pct))
                    st.session_state.nominal_sizes = st.session_state.nominal_sizes.append({'%F(d)':pct,'Tamaño (µm)':d_val}, ignore_index=True)
                    st.success(f"Grabado: {pct}% -> {d_val:.3f} µm")
                except Exception as e:
                    st.error("No se puede interpolar: " + str(e))

        if st.button("CALCULAR PORCENTAJES PASANTES"):
            d_input = st.number_input("Ingrese d (µm) =", min_value=0.0, value=100.0, key='calc_pct_d')
            if st.button("GRABAR d->%F(d)", key='grabar_pct'):
                x = results['Tamaño promedio (µm)'].replace(0, np.nan)
                y = results['%F(d)']
                mask = ~np.isnan(x) & ~np.isnan(y)
                try:
                    f = interp1d(x[mask], y[mask], fill_value="extrapolate")
                    pctv = float(f(d_input))
                    st.session_state.nominal_sizes = st.session_state.nominal_sizes.append({'%F(d)':pctv,'Tamaño (µm)':d_input}, ignore_index=True)
                    st.success(f"Grabado: {d_input} µm -> {pctv:.3f}%")
                except Exception as e:
                    st.error("No se puede interpolar: " + str(e))

        st.markdown("**Tabla de tamaños nominales grabados**")
        st.dataframe(st.session_state.nominal_sizes)

    # ---------- Sección 3: Estadísticos según Folk & Ward ----------
    st.header("ESTADÍSTICOS SEGÚN FOLK Y WARD")
    st.markdown("Se usan tamaños nominales grabados en la sección anterior. Si no existen los tamaños requeridos, se mostrará un warning.")
    required_perc = [5,16,25,50,75,84]
    nom_df = st.session_state.nominal_sizes.copy()
    if nom_df.empty:
        st.warning("Calcule los tamaños nominales requeridos en la sección anterior para mostrar los estadísticos.")
    else:
        x = results['Tamaño promedio (µm)'].replace(0, np.nan)
        y = results['%F(d)']
        mask = ~np.isnan(x) & ~np.isnan(y)
        try:
            inv = interp1d(y[mask], x[mask], fill_value="extrapolate")
            ds = {}
            for p in required_perc:
                ds[p] = float(inv(p))
            phi = {p: -np.log2(ds[p]) for p in ds}
            M = (phi[16] + phi[50] + phi[84]) / 3.0
            Md = phi[50]
            sigma = (phi[84] - phi[16]) / 4.0
            Sk = (phi[16] + phi[84] - 2*phi[50]) / (2*(phi[84]-phi[16])) if (phi[84]-phi[16])!=0 else 0.0
            K = np.nan
            folk_tbl = pd.DataFrame({
                'Nombre':['M (media, phi)','Md (mediana, phi)','σ (varianza/dispersion, phi)','Sk (asimetría)','K (curtosis)'],
                'Valor':[M, Md, sigma, Sk, K]
            })
            st.table(folk_tbl)
            if sigma > 1:
                disp_comment = "La muestra presenta una amplia dispersión de tamaños."
            else:
                disp_comment = "La muestra presenta una dispersión reducida de tamaños."
            skew_comment = "La distribución es sesgada hacia partículas gruesas." if Sk < 0 else "La distribución es sesgada hacia partículas finas."
            kurt_comment = "Curtosis alta (leptocúrtica) o baja según valor K."
            st.markdown(f"**Interpretación:** {disp_comment} {skew_comment} {kurt_comment}")
        except Exception as e:
            st.error("No se pueden calcular Folk & Ward: " + str(e))

    if st.button("SIGUIENTE"):
        st.session_state.page = 5
        st.rerun()

# ---------- PÁGINA 5: Selección del modelo ----------
def page_5():
    st.title("SELECCIÓN DEL MODELO")
    st.markdown("Se introducirán y ajustarán los modelos: GGS, RRSB y Doble Weibull. Se estiman parámetros minimizando la función objetivo (SSE).")
    results = st.session_state.results_table.copy()
    if results.empty:
        st.error("No hay datos procesados.")
        if st.button("Regresar"):
            st.session_state.page = 4
            st.rerun()
        return
    d = results['Tamaño promedio (µm)'].astype(float).values
    y_exp = results['%F(d)'].astype(float).values
    mask = (d > 0) & (~np.isnan(y_exp))
    d_fit = d[mask]
    y_fit = y_exp[mask]

    if len(d_fit) < 3:
        st.warning("No hay suficientes puntos válidos para ajuste (se requieren al menos 3).")
    else:
        if st.button("AJUSTAR"):
            x0 = [1.0, np.median(d_fit)]
            def f_ggs(params):
                m, Dm = params
                return np.sum((y_fit - GGS_model(d_fit, m, Dm))**2)
            res1 = minimize(f_ggs, x0, bounds=[(0.01,10),(1e-6, max(d_fit)*10)])
            FO_ggs = res1.fun
            ggs_params = res1.x

            x0 = [1.0, np.median(d_fit)]
            def f_rrsb(params):
                m, l = params
                return np.sum((y_fit - RRSB_model(d_fit, m, l))**2)
            res2 = minimize(f_rrsb, x0, bounds=[(0.01,10),(1e-6, max(d_fit)*10)])
            FO_rrsb = res2.fun
            rrsb_params = res2.x

            try:
                x_all = results['Tamaño promedio (µm)'].replace(0, np.nan)
                y_all = results['%F(d)']
                mask2 = ~np.isnan(x_all) & ~np.isnan(y_all)
                inv = interp1d(y_all[mask2], x_all[mask2], fill_value="extrapolate")
                init_d80 = float(inv(80.0))
            except Exception:
                init_d80 = np.median(d_fit)
            x0 = [0.5, 1.0, 1.0, init_d80]
            def f_double(params):
                alpha, k1, k2, d80 = params
                return np.sum((y_fit - 
double_weibull(d_fit, alpha, k1, k2, d80))**2)
            bounds_dw = [(0.0,1.0),(0.01,10.0),(0.01,10.0),(1e-3,max(d_fit)*10)]
            res3 = minimize(f_double, x0, bounds=bounds_dw)
            FO_dw = res3.fun
            dw_params = res3.x

            st.session_state.models_fit = {
                'GGS': {'FO':FO_ggs, 'params':ggs_params},
                'RRSB': {'FO':FO_rrsb, 'params':rrsb_params},
                'DobleWeibull': {'FO':FO_dw, 'params':dw_params}
            }
            st.success("Ajustes realizados.")
            st.rerun()

        if st.session_state.models_fit:
            fits = st.session_state.models_fit
            st.subheader("Comparación de funciones objetivo (FO)")
            fo_tbl = pd.DataFrame([
                {'Modelo':'GGS','F.O.':fits['GGS']['FO']},
                {'Modelo':'RRSB','F.O.':fits['RRSB']['FO']},
                {'Modelo':'Doble Weibull','F.O.':fits['DoubleWeibull']['FO']}
            ])
            st.table(fo_tbl)
            best = min(fits.items(), key=lambda x: x[1]['FO'])
            best_model_name = best[0]
            st.markdown(f"**El mejor modelo es {best_model_name} con F.O.= {best[1]['FO']:.3f}**")
            if best[1]['FO'] > 1e6:
                st.warning("Ningún modelo representa bien los datos experimentales (F.O. muy grande).")
            st.subheader("Parámetros estimados del mejor modelo")
            st.json(best[1]['params'])
            
            fig, ax = plt.subplots(figsize=(7,4))
            aplicar_estilo_granulometrico(ax)
            ax.plot(d_fit, y_fit, 'o', label='Experimental', mfc='white', mec='#004488', color='#004488', ms=8)
            dd = np.linspace(np.min(d_fit), np.max(d_fit), 200)
            if best_model_name == 'GGS':
                m, Dm = best[1]['params']
                ax.plot(dd, GGS_model(dd, m, Dm), '-', label=f'GGS (m={m:.3f}, Dm={Dm:.3f})', color='#004488', lw=2)
            elif best_model_name == 'RRSB':
                m, l = best[1]['params']
                ax.plot(dd, RRSB_model(dd, m, l), '-', label=f'RRSB (m={m:.3f}, l={l:.3f})', color='#004488', lw=2)
            else:
                alpha, k1, k2, d80 = best[1]['params']
                ax.plot(dd, double_weibull(dd, alpha, k1, k2, d80), '-', label=f'DW (alpha={alpha:.3f},k1={k1:.3f},k2={k2:.3f},d80={d80:.3f})', color='#004488', lw=2)
            ax.set_xlabel("Tamaño (µm)")
            ax.set_ylabel("%F(d)")
            ax.legend()
            ax.grid(True)
            st.pyplot(fig)

    if st.button("SIGUIENTE"):
        st.session_state.page = 6
        st.rerun()

# ---------- PÁGINA 6: Exportación ----------
def page_6():
    st.title("EXPORTACIÓN DE DATOS")
    st.markdown("Descargar todas las tablas (sin gráficos) en un archivo Excel o guardar el análisis y generar un QR para compartir.")

    # Prepare Excel in memory
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        st.session_state.input_table.to_excel(writer, sheet_name='Entrada', index=False)
        st.session_state.results_table.to_excel(writer, sheet_name='Resultados', index=False)
        st.session_state.nominal_sizes.to_excel(writer, sheet_name='TamañosNominales', index=False)
        models = st.session_state.models_fit
        if models:
            models_df = pd.DataFrame([
                {'Modelo':k, 'FO':v['FO'], 'Parametros':str(v['params'])} for k,v in models.items()
            ])
            models_df.to_excel(writer, sheet_name='Modelos', index=False)
    data = output.getvalue()

    b64 = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{b64}" download="analisis_granulometrico.xlsx">Descargar Excel (analisis_granulometrico.xlsx)</a>'
    st.markdown(href, unsafe_allow_html=True)

    if st.button("GUARDAR (generar QR de descarga)"):
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx')
        with open(tmp.name, 'wb') as f:
            f.write(data)
        link = f"file://{tmp.name}"
        download_url = f"data:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;base64,{b64}"
        qr = qrcode.make(download_url)
        st.image(qr)
        st.markdown("Escanea el QR para descargar el archivo (si tu lector soporta data URLs).")
        st.success(f"Archivo guardado temporalmente en: {tmp.name}")

    if st.button("VOLVER AL INICIO"):
        st.session_state.page = 1
        st.rerun()

# ---------- Page router ----------
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
    else:
        st.session_state.page = 1
        page_1()

if __name__ == "__main__":
    main()












