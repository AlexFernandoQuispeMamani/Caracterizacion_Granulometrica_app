import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import xlsxwriter
import math
from scipy.optimize import curve_fit

# Definir la serie de mallas Tyler Standard
tyler_series = {
    '"1.05"': 26500, '"0.883"': 22400, '"0.742"': 19000, '"0.624"': 16000,
    '"0.525"': 13200, '"0.441"': 11200, '"0.371"': 9500, '2 ½”': 8000,
    '3”': 6700, '3 ½”': 5600, '4#': 4760, '5#': 4000, '6#': 3350,
    '7#': 2800, '8#': 2360, '9#': 2000, '10#': 1700, '12#': 1400,
    '14#': 1180, '16#': 1000, '20#': 850, '24#': 710, '28#': 600,
    '32#': 500, '35#': 425, '42#': 355, '48#': 300, '60#': 250,
    '65#': 212, '80#': 180, '100#': 150, '115#': 125, '150#': 106,
    '170#': 90, '200#': 75, '250#': 63, '270#': 53, '320#': 45,
    '400#': 38
}

def ggs_model(x, k, m):
    """Modelo de Gates-Gaudin-Schuhmann (GGS)"""
    return 100 * (x / k)**m

def rrsb_model(x, d_prime, n):
    """Modelo de Rosin-Rammler-Sperling-Bennett (RRSB)"""
    return 100 * (1 - np.exp(-(x / d_prime)**n))

def weibull_model(x, alpha, beta, gamma):
    """Modelo de Weibull"""
    return 100 * (1 - np.exp(-((x - gamma) / alpha)**beta))

def create_and_download_excel(dfs, person_info, plot_data_bytes):
    """Crea un archivo Excel en memoria con múltiples hojas."""
    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        
        # Guardar la información personal en la primera hoja
        person_df = pd.DataFrame([person_info]).T.rename(columns={0: 'Valor'})
        person_df.index.name = 'Dato'
        person_df.to_excel(writer, sheet_name='Información General')
        
        # Guardar cada DataFrame en una hoja separada
        for sheet_name, df in dfs.items():
            df.to_excel(writer, sheet_name=sheet_name)
            
        # Insertar los gráficos
        worksheet = writer.sheets['Análisis Granulométrico']
        if plot_data_bytes:
            worksheet.insert_image('G2', 'imagen.png', {'image_data': plot_data_bytes})
    
    return output.getvalue()

st.set_page_config(layout="wide")

if 'page' not in st.session_state:
    st.session_state.page = "home"
    st.session_state.user_data = {}
    st.session_state.exp_data = {}
    st.session_state.results = {}

def set_page(page_name):
    st.session_state.page = page_name

def validate_home_data():
    if (st.session_state.user_data.get('usuario') and
        st.session_state.user_data.get('correo') and
        st.session_state.user_data.get('procedencia') and
        st.session_state.user_data.get('codigo') and
        st.session_state.user_data.get('fecha')):
        return True
    return False

def page_1_home():
    st.title("CARACTERIZACIÓN GRANULOMÉTRICA")
    st.markdown("""
    Desarrollado por **Alex Fernando Quispe Mamani**, ingeniero metalúrgico graduado de la UTO con experiencia en desarrollo y simulación de procesos metalúrgicos y escritor de textos de pregrado.

    Este programa realiza un análisis granulométrico completo de una muestra de mineral, utilizando datos experimentales de tamizado. Permite al usuario introducir los datos de diferentes maneras y calcula la distribución de tamaño de partícula, los estadísticos de Folk y Ward, y ajusta los datos a los modelos de Gates-Gaudin-Schuhmann (GGS), Rosin-Rammler-Sperling-Bennett (RRSB) y Weibull.

    ---

    ### INFORMACIÓN GENERAL
    """)
    with st.form("personal_info_form"):
        st.session_state.user_data['usuario'] = st.text_input("USUARIO (su nombre)")
        st.session_state.user_data['correo'] = st.text_input("CORREO ELECTRÓNICO (gmail)")
        st.session_state.user_data['procedencia'] = st.text_input("PROCEDENCIA DE LA MUESTRA")
        st.session_state.user_data['codigo'] = st.text_input("CÓDIGO DE LA MUESTRA")
        st.session_state.user_data['fecha'] = st.date_input("FECHA DE MUESTREO")
        
        submitted = st.form_submit_button("INICIO")
        if submitted:
            if validate_home_data():
                set_page("experimental_data")
            else:
                st.warning("Por favor, complete todos los campos para continuar.")

def page_2_experimental_data():
    st.title("DATOS EXPERIMENTALES")
    st.write("Por favor, inserte los tamaños y pesos retenidos sobre cada tamiz para efectuar el análisis granulométrico.")
    st.write("Seleccione **SELECCIONAR MALLAS** si tiene el número de malla de los tamices utilizados o **INSERTAR MANUALMENTE** para insertar los tamaños en micrones de manera personalizada.")
    
    st.session_state.exp_data['tipo_entrada'] = st.radio("Seleccione el tipo de entrada de datos:", ["SELECCIONAR MALLAS", "INSERTAR MANUALMENTE"])

    if st.session_state.exp_data['tipo_entrada'] == "SELECCIONAR MALLAS":
        st.session_state.exp_data['peso_total'] = st.number_input("Peso total (g):", min_value=0.0)
        selected_sieves = st.multiselect("Seleccionar Mallas:", options=list(tyler_series.keys()), format_func=lambda x: f"{x} ({tyler_series[x]} µm)")
        
        if selected_sieves:
            selected_sieves_sorted = sorted(selected_sieves, key=lambda x: tyler_series[x], reverse=True)
            sizes = [tyler_series[s] for s in selected_sieves_sorted]
            weights = [0.0] * len(selected_sieves_sorted)
            data = {'Nº de malla': selected_sieves_sorted, 'Tamaño (μm)': sizes, 'Peso (g)': weights}
            df = pd.DataFrame(data)
            st.session_state.exp_data['df'] = st.data_editor(df, num_rows="dynamic", hide_index=True)
    
    else: # INSERTAR MANUALMENTE
        st.session_state.exp_data['peso_total'] = st.number_input("Peso total (g):", min_value=0.0)
        num_rows = st.number_input("Número de datos a insertar (3-25):", min_value=3, max_value=25, value=3)
        df_manual = pd.DataFrame({'Tamaño (μm)': [0.0] * num_rows, 'Peso (g)': [0.0] * num_rows})
        st.session_state.exp_data['df'] = st.data_editor(df_manual, num_rows="dynamic", hide_index=True)
    
    st.info("Recomendación: Asegúrese de que los tamaños de los tamices estén ordenados de mayor a menor.")
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            set_page("home")
    with col2:
        if st.button("EJECUTAR"):
            if 'df' in st.session_state.exp_data and not st.session_state.exp_data['df'].empty and st.session_state.exp_data['peso_total'] > 0:
                set_page("analysis_results")
            else:
                st.warning("Por favor, ingrese datos válidos en la tabla y un peso total para continuar.")
                
def page_3_analysis_results():
    st.title("ANÁLISIS GRANULOMÉTRICO")

    df_data = st.session_state.exp_data['df'].copy()
    
    # Cálculos
    df_data['Peso retenido acumulado (g)'] = df_data['Peso (g)'].cumsum()
    df_data['% Retenido'] = (df_data['Peso (g)'] / df_data['Peso (g)'].sum()) * 100
    df_data['% Retenido acumulado'] = df_data['% Retenido'].cumsum()
    df_data['% Pasa'] = 100 - df_data['% Retenido acumulado']

    st.session_state.results['df_granulometrico'] = df_data
    
    st.write("### Tabla de Análisis Granulométrico")
    
    if st.session_state.exp_data['tipo_entrada'] == "SELECCIONAR MALLAS":
        # Lógica para la columna 'Nº de malla' en forma de intervalo
        mallas = st.session_state.exp_data['df']['Nº de malla'].tolist()
        malla_intervalo = []
        if len(mallas) > 0:
            malla_intervalo.append(f"{mallas[0]}") # Primera malla
            for i in range(len(mallas) - 1):
                malla_intervalo.append(f"-{mallas[i]}+{mallas[i+1]}")
            malla_intervalo.append(f"-{mallas[-1]}")
        
        df_intervalo = pd.DataFrame({'Nº de malla': malla_intervalo, '% Pasa': df_data['% Pasa'].tolist() + [0.0]})
        df_intervalo.loc[len(df_intervalo)-1, '% Pasa'] = 0.0 # Asegurar que el último valor es 0
        df_intervalo = df_intervalo.set_index('Nº de malla')
        st.table(df_intervalo)

    st.table(df_data)

    st.write("---")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            set_page("experimental_data")
    with col2:
        if st.button("SIGUIENTE"):
            set_page("statistical_analysis")

def page_4_statistical_analysis():
    st.title("ANÁLISIS DE DATOS")

    df_gran = st.session_state.results['df_granulometrico']
    
    # Se obtienen los valores para los estadísticos
    def get_d(df, p):
        df_sorted = df.sort_values(by='Tamaño (μm)', ascending=False)
        passing_sizes = df_sorted['% Pasa'].tolist()
        sizes = df_sorted['Tamaño (μm)'].tolist()
        
        if p > max(passing_sizes) or p < min(passing_sizes):
            # Extrapolación log-log
            # Para puntos que pasan por encima del máximo
            if p > max(passing_sizes) and len(sizes) > 1:
                x1, y1 = np.log10(sizes[0]), np.log10(passing_sizes[0])
                x2, y2 = np.log10(sizes[1]), np.log10(passing_sizes[1])
                x_interp = np.log10(p)
                y_interp = y1 + (x_interp - y1) * (x2 - x1) / (y2 - y1) # Falla la formula?
                # Interpolación con regresión lineal
                log_pasa = np.log10(df_sorted['% Pasa'])
                log_size = np.log10(df_sorted['Tamaño (μm)'])
                m, c = np.polyfit(log_pasa, log_size, 1)
                log_d = m * np.log10(p) + c
                return 10**log_d
                
            # Para puntos que pasan por debajo del mínimo
            elif p < min(passing_sizes) and len(sizes) > 1:
                x1, y1 = np.log10(sizes[-2]), np.log10(passing_sizes[-2])
                x2, y2 = np.log10(sizes[-1]), np.log10(passing_sizes[-1])
                x_interp = np.log10(p)
                y_interp = y1 + (x_interp - y1) * (x2 - x1) / (y2 - y1) # Falla la formula?
                # Interpolación con regresión lineal
                log_pasa = np.log10(df_sorted['% Pasa'])
                log_size = np.log10(df_sorted['Tamaño (μm)'])
                m, c = np.polyfit(log_pasa, log_size, 1)
                log_d = m * np.log10(p) + c
                return 10**log_d
            else:
                return np.nan # No se puede extrapolar con menos de 2 puntos
        
        return np.interp(p, passing_sizes, sizes)

    # Cálculo de los estadísticos
    try:
        d16 = get_d(df_gran, 16)
        d50 = get_d(df_gran, 50)
        d84 = get_d(df_gran, 84)

        media = (d16 + d50 + d84) / 3
        desviacion_estandar = ((d84 - d16) / 4) + ((d95 - d5) / 6.6)
        
        estadisticos = {
            "d16": d16, "d50": d50, "d84": d84,
            "Media (ϕ)": media,
            "Desviación Estándar (ϕ)": desviacion_estandar,
        }

        st.session_state.results['estadisticos'] = pd.DataFrame(estadisticos.items(), columns=['Estadístico', 'Valor']).set_index('Estadístico')
        st.write("### Estadísticos de la Muestra")
        st.table(st.session_state.results['estadisticos'])
    except:
        st.warning("No se pudieron calcular los estadísticos debido a la falta de datos.")
    
    st.write("---")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            set_page("analysis_results")
    with col2:
        if st.button("SIGUIENTE"):
            set_page("model_selection")

def page_5_model_selection():
    st.title("SELECCIÓN DEL MODELO")

    df_gran = st.session_state.results['df_granulometrico']
    x_data = df_gran['Tamaño (μm)'].values
    y_data = df_gran['% Pasa'].values
    
    try:
        # Ajuste de modelos
        params_ggs, _ = curve_fit(ggs_model, x_data, y_data)
        params_rrsb, _ = curve_fit(rrsb_model, x_data, y_data)
        
        # Generar puntos para las curvas de los modelos
        x_fit = np.logspace(np.log10(min(x_data)), np.log10(max(x_data)), 100)
        y_ggs = ggs_model(x_fit, *params_ggs)
        y_rrsb = rrsb_model(x_fit, *params_rrsb)

        # Calcular Función Objetivo (F.O.)
        def calculate_fo(model_func, params, x, y):
            y_pred = model_func(x, *params)
            return np.sum((y_pred - y)**2)

        fo_ggs = calculate_fo(ggs_model, params_ggs, x_data, y_data)
        fo_rrsb = calculate_fo(rrsb_model, params_rrsb, x_data, y_data)
        
        st.session_state.results['fo_values'] = {
            "GGS": fo_ggs,
            "RRSB": fo_rrsb
        }
        
        # Gráficos
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(x_data, y_data, 'o', label='Datos Experimentales')
        ax.plot(x_fit, y_ggs, '-', label='Modelo GGS')
        ax.plot(x_fit, y_rrsb, '-', label='Modelo RRSB')
        ax.set_xscale('log')
        ax.set_ylim(0, 100)
        ax.set_xlabel('Tamaño (μm)')
        ax.set_ylabel('% Pasa')
        ax.set_title('Curvas de Distribución Granulométrica')
        ax.grid(True, which="both", ls="--")
        ax.legend()
        st.pyplot(fig)
        
        # Tabla comparativa y mejor modelo
        st.write("### Valores de la Función Objetivo (F.O.)")
        df_fo = pd.DataFrame(st.session_state.results['fo_values'].items(), columns=['Modelo', 'F.O.']).set_index('Modelo')
        st.table(df_fo)
        
        best_model = min(st.session_state.results['fo_values'], key=st.session_state.results['fo_values'].get)
        st.write(f"El mejor modelo que se ajusta a los datos es: **{best_model}** con un valor de F.O. de {st.session_state.results['fo_values'][best_model]:.2f}.")
    except Exception as e:
        st.error(f"Error al ajustar los modelos: {e}. Asegúrese de que sus datos son válidos.")

    st.write("---")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("ANTERIOR"):
            set_page("statistical_analysis")
    with col2:
        if st.button("SIGUIENTE"):
            set_page("export_data")

def page_6_export_data():
    st.title("EXPORTACIÓN DE DATOS")
    
    st.write("Haga clic en el botón 'DESCARGA' para descargar un archivo Excel con todas las tablas de su análisis.")

    # Convertir las tablas en dataframes
    all_dfs = {
        "Información General": pd.DataFrame([st.session_state.user_data]).T.rename(columns={0: 'Valor'}),
        "Datos Experimentales": st.session_state.exp_data['df'],
        "Análisis Granulométrico": st.session_state.results['df_granulometrico'],
        "Estadísticos": st.session_state.results.get('estadisticos', pd.DataFrame())
    }
    
    # Crear la imagen del gráfico para exportar
    fig, ax = plt.subplots(figsize=(10, 6))
    if 'df_granulometrico' in st.session_state.results:
        df_gran = st.session_state.results['df_granulometrico']
        ax.plot(df_gran['Tamaño (μm)'], df_gran['% Pasa'], 'o', label='Datos Experimentales')
        ax.set_xscale('log')
        ax.set_ylim(0, 100)
        ax.set_xlabel('Tamaño (μm)')
        ax.set_ylabel('% Pasa')
        ax.set_title('Curva Granulométrica')
        ax.grid(True, which="both", ls="--")
        ax.legend()
    
    img_buffer = BytesIO()
    fig.savefig(img_buffer, format='png')
    img_buffer.seek(0)
    
    excel_data = create_and_download_excel(all_dfs, st.session_state.user_data, img_buffer)
    
    st.download_button(
        label="DESCARGA",
        data=excel_data,
        file_name='Analisis_Granulometrico.xlsx',
        mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )
    
    st.write("---")
    if st.button("ANTERIOR"):
        set_page("model_selection")
        
# Lógica para mostrar la página actual
if st.session_state.page == "home":
    page_1_home()
elif st.session_state.page == "experimental_data":
    page_2_experimental_data()
elif st.session_state.page == "analysis_results":
    page_3_analysis_results()
elif st.session_state.page == "statistical_analysis":
    page_4_statistical_analysis()
elif st.session_state.page == "model_selection":
    page_5_model_selection()
elif st.session_state.page == "export_data":
    page_6_export_data()




