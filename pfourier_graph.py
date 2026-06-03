import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import config as cf  # Importamos tu archivo de configuración

# ==========================================
# CONFIGURACIÓN: Rutas y Datasets
# ==========================================
ARCHIVO_HDF5 = "wave_data.hdf5"
DATASET_DFT = "/Probes/Probe1"  
NOMBRE_GRAFICO_SALIDA = "respuesta_frecuencia_dft.pdf" 

# --- CONTROL MANUAL DE UNIDADES DEL EJE X ---
# Configura aquí el factor de escala y la etiqueta que desees:
DIVISOR_X = 1e3      # Usa: 1 para Hz, 1e6 para MHz, 1e9 para GHz, etc.
UNIDAD_X = "KHz"     # El texto que acompañará al eje X: "Hz", "MHz", "GHz"...

# Constante de la luz en el vacío (m/s)
C_LIGHT = 299792458.0 
# ==========================================

# ==========================================
# ESTILO DE PUBLICACIÓN (Paper Style)
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

def graficar_dft_fisica(ruta_archivo):
    """Lee la DFT y grafica la respuesta en frecuencia con escala e unidades personalizables."""
    try:
        with h5py.File(ruta_archivo, "r") as f:
            datos_dft = f[DATASET_DFT][:]

        # ---------------------------------------------------------
        # FÍSICA Y MATEMÁTICAS: Cálculo de resoluciones
        # ---------------------------------------------------------
        lambda_min = C_LIGHT / cf.MAX_FREQ
        delta_x = lambda_min / 20.0
        delta_t = cf.COURANT * (delta_x / C_LIGHT)

        t_total = cf.TOTAL_TIME * delta_t
        f_min = 1.0 / t_total
        delta_f = f_min

        # Vector de frecuencias original (en Hz)
        frecuencias_hz = np.arange(len(datos_dft)) * delta_f
        
        # Aplicamos tu división manual para la visualización
        frecuencias_escaladas = frecuencias_hz / DIVISOR_X
        
        # Procesamiento del DFT (solo magnitud)
        magnitud = np.abs(datos_dft)

        # ---------------------------------------------------------
        # GENERACIÓN DEL GRÁFICO
        # ---------------------------------------------------------
        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)

        # Trazado de la señal con los datos escalados
        ax.plot(
            frecuencias_escaladas, 
            magnitud, 
            color="black", 
            linewidth=1.5,
            zorder=3
        )

        # Configuración dinámica de los ejes
        ax.set_xlabel(f"Frecuencia ({UNIDAD_X})", fontsize=13)
        ax.set_ylabel("Espectro normalizado", fontsize=13)
        
        # El límite también se ajusta automáticamente con tu divisor
        ax.set_xlim(0, cf.MAX_FREQ / DIVISOR_X)
        ax.set_ylim(bottom=0)

        # Forzamos a que no meta multiplicadores raros por su cuenta
        ax.ticklabel_format(style='plain', axis='x', useOffset=False)

        ax.grid(False) 
        plt.tight_layout()
        
        # ---------------------------------------------------------
        # GUARDAR LA IMAGEN
        # ---------------------------------------------------------
        plt.savefig(f"./fotogramas/{NOMBRE_GRAFICO_SALIDA}", dpi=600, bbox_inches='tight')
        print(f"[INFO] Imagen guardada con éxito como: {NOMBRE_GRAFICO_SALIDA}")

        plt.show()

        # Resumen en consola para verificación rápida
        print("\n--- Parámetros Físicos Calculados ---")
        print(f"Malla espacial (delta_x) : {delta_x:.4e} m")
        print(f"Paso de tiempo (delta_t) : {delta_t:.4e} s")
        print(f"Tiempo de simulación     : {t_total:.4e} s")
        print(f"Resolución (delta_f)     : {delta_f:.2f} Hz")
        print(f"Frecuencia máxima (lim)  : {cf.MAX_FREQ:.2e} Hz")
        print("-------------------------------------")

    except KeyError:
        print(f"\n[ERROR]: No se encontró el dataset '{DATASET_DFT}' en el archivo HDF5.")
    except AttributeError as e:
        print(f"\n[ERROR]: Falta una variable en tu archivo config.py. Detalles: {e}")
    except FileNotFoundError:
        print(f"\n[ERROR]: El archivo '{ruta_archivo}' no existe en este directorio.")

if __name__ == "__main__":
    graficar_dft_fisica(ARCHIVO_HDF5)