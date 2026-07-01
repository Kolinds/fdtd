import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import config as cf

# ==========================================
# CONFIGURACIÓN
# ==========================================
# Definimos los 6 archivos y sus etiquetas para la leyenda
ARCHIVOS_HDF5 = ["reflec_trans1.hdf5", "reflec_trans4.hdf5", "reflec_trans5.hdf5", "reflec_trans6.hdf5"]

ETIQUETAS = ["ADE", "Implicit $S_i = S_c$", "Implicit $S_i = 2S_c$", "Implicit $S_i = 3S_c$"]

# Paleta corregida sin espacios para evitar el ValueError
COLORES = ["black", "firebrick", "royalblue", "forestgreen", "darkorchid", "darkslategray"]

DATASET_TRANSMISSION = "reflection"

# --- INTERRUPTOR DE ESCALA ---
ESCALA_DB = False  # True para escala logarítmica (dB), False para escala lineal
# -----------------------------

# --- MARCADOR DE FRECUENCIA DE PLASMA ---
MOSTRAR_FREQ_PLASMA = True  # Cambia a False para ocultar la línea
# El valor ahora se calcula automáticamente usando cf.PLASMA_WAVELENGTH_STEPS
# ----------------------------------------

# Nombre base del gráfico (se le añadirá _db o _lineal automáticamente)
NOMBRE_GRAFICO_BASE = "transmitancia_multiples_casos"

# Escala del eje X
DIVISOR_X = 1e6      # 1=Hz, 1e3=kHz, 1e6=MHz, 1e9=GHz
UNIDAD_X = "MHz"

# Velocidad de la luz (m/s)
C_LIGHT = 299792458.0
# ==========================================

# ==========================================
# ESTILO DE PUBLICACIÓN
# ==========================================
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.major.width'] = 1.2
mpl.rcParams['ytick.major.width'] = 1.2
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
# ==========================================


def graficar_transmitancias_multiples(lista_archivos, lista_etiquetas):
    """
    Lee el coeficiente de transmisión complejo de múltiples archivos HDF5,
    grafica los resultados y calcula automáticamente fp a partir de los pasos de malla.
    """
    # Inicializamos la gráfica fuera del bucle para que todo vaya al mismo lienzo
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)
    
    longitud_puntos = 0
    delta_f_control = 0
    delta_x_control = 0
    delta_t_control = 0
    t_total_control = 0

    # Iteramos sobre cada archivo
    for i, ruta_archivo in enumerate(lista_archivos):
        try:
            # -------------------------------------------------
            # CARGA DE DATOS
            # -------------------------------------------------
            with h5py.File(ruta_archivo, "r") as f:
                t_complex = f[DATASET_TRANSMISSION][:]

            # -------------------------------------------------
            # PARÁMETROS FDTD
            # -------------------------------------------------
            lambda_min = C_LIGHT / cf.MAX_FREQ
            delta_x = lambda_min / 20.0
            delta_t = cf.COURANT * (delta_x / C_LIGHT)
            t_total = cf.TOTAL_TIME * delta_t
            delta_f = 1.0 / t_total

            # Guardamos variables del último archivo para el reporte final de consola
            longitud_puntos = len(t_complex)
            delta_f_control = delta_f
            delta_x_control = delta_x
            delta_t_control = delta_t
            t_total_control = t_total

            # -------------------------------------------------
            # EJE DE FRECUENCIAS Y TRANSMITANCIA
            # -------------------------------------------------
            frecuencias_hz = np.arange(len(t_complex)) * delta_f
            frecuencias = frecuencias_hz / DIVISOR_X
            
            # Magnitud lineal absoluta del coeficiente complejo
            magnitud_lineal = np.abs(t_complex)
            
            # Selección de datos a graficar según el interruptor
            if ESCALA_DB:
                datos_y = 20 * np.log10(magnitud_lineal + 1e-12)
            else:
                datos_y = magnitud_lineal

            # -------------------------------------------------
            # AÑADIR LÍNEA A LA GRÁFICA
            # -------------------------------------------------
            ax.plot(
                frecuencias,
                datos_y,
                color=COLORES[i],
                label=lista_etiquetas[i],
                linewidth=1.5,
                zorder=3
            )
            print(f"[OK] Datos cargados para: {ruta_archivo}")

        except KeyError:
            print(f"[ERROR] No se encontró el dataset '{DATASET_TRANSMISSION}' en '{ruta_archivo}'.")
            return
        except FileNotFoundError:
            print(f"[ERROR] El archivo '{ruta_archivo}' no existe.")
            return
        except AttributeError as e:
            print(f"[ERROR] Falta alguna variable en config.py:\n{e}")
            return

    # ==========================================
    # CONFIGURACIÓN FINAL DEL LIENZO (Post-bucle)
    # ==========================================
    
    # --- CÁLCULO DINÁMICO Y LÍNEA DE FRECUENCIA DE PLASMA ---
    if MOSTRAR_FREQ_PLASMA:
        try:
            # lambda_p = pasos_de_malla * delta_x
            lambda_plasma = cf.PLASMA_WAVELENGTH_STEPS * delta_x_control
            # f_p = c / lambda_p
            freq_plasma_hz = C_LIGHT / lambda_plasma
            # Convertimos a la unidad del eje (MHz o GHz)
            freq_plasma_grafico = freq_plasma_hz / DIVISOR_X

            ax.axvline(
                x=freq_plasma_grafico, 
                color='darkorange', 
                linewidth=1.5, 
                linestyle='--', 
                label=rf'$f_p$',
                zorder=4
            )
        except AttributeError:
            print("[ERROR] No se encontró 'PLASMA_WAVELENGTH_STEPS' en config.py")
    # -------------------------------------

    ax.set_xlabel(f"Frecuencia ({UNIDAD_X})", fontsize=13)
    
    # Ajuste de límites X
    ax.set_xlim(0.1, 0.7)

    # Ajuste condicional del eje Y según el interruptor
    if ESCALA_DB:
        ax.set_ylabel("Coef. transmisión (dB)", fontsize=13)
        ax.set_ylim(-80, 2.0)
        nombre_salida = f"{NOMBRE_GRAFICO_BASE}_db.pdf"
    else:
        ax.set_ylabel("Coef. reflexión (u.a.)", fontsize=13)
        ax.set_ylim(-0.05, 1.1)
        nombre_salida = f"{NOMBRE_GRAFICO_BASE}_lineal.pdf"

    ax.ticklabel_format(style='plain', axis='x', useOffset=False)
    
    # Añadimos la leyenda (ajustada automáticamente)
    ax.legend(loc="best", frameon=True, edgecolor="none", fontsize=10)
    
    # Cuadrícula tenue
    ax.grid(True, linestyle=":", alpha=0.5, color="gray")
    plt.tight_layout()

    # -------------------------------------------------
    # GUARDAR Y MOSTRAR FIGURA
    # -------------------------------------------------
    plt.savefig(f"./fotogramas/{nombre_salida}", dpi=600, bbox_inches='tight')
    print(f"\n[INFO] Imagen multienfoque guardada con éxito como: {nombre_salida}")
    plt.show()

    # -------------------------------------------------
    # INFORMACIÓN DE CONTROL (Común a las simulaciones)
    # -------------------------------------------------
    print("\n--- Parámetros FDTD Compartidos ---")
    print(f"delta_x           = {delta_x_control:.4e} m")
    print(f"delta_t           = {delta_t_control:.4e} s")
    print(f"Tiempo total      = {t_total_control:.4e} s")
    print(f"delta_f           = {delta_f_control:.4e} Hz")
    if MOSTRAR_FREQ_PLASMA:
        print(f"Frecuencia plasma = {freq_plasma_hz:.4e} Hz ({freq_plasma_grafico:.2f} {UNIDAD_X})")
    print(f"Frecuencia máxima = {cf.MAX_FREQ:.4e} Hz")
    print(f"Número de puntos  = {longitud_puntos}")
    print("-----------------------------------")


if __name__ == "__main__":
    graficar_transmitancias_multiples(ARCHIVOS_HDF5, ETIQUETAS)