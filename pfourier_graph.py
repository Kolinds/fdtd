import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import config as cf  # Importamos tu archivo de configuración

# ==========================================
# CONFIGURACIÓN: Rutas y Datasets
# ==========================================
ARCHIVO_HDF5 = "wave_data.hdf5"
DATASET_DFT = "/Probes/Reflected"   

NOMBRE_GRAFICO_SALIDA = "respuesta_frecuencia_dft.pdf" 

# --- CONTROL DE LAS LÍNEAS VERTICALES ---
MOSTRAR_LINEA_PLASMA = True  # Activa/Desactiva la línea de frecuencia de plasma
MOSTRAR_LINEA_PEAK = True    # ¡NUEVO! Activa/Desactiva la línea de frecuencia pico

# --- CONTROL MANUAL DE UNIDADES DEL EJE X ---
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
        frecuencias_escaladas = frecuencias_hz / DIVISOR_X
        
        # Cálculo de la frecuencia de plasma
        lambda_plasma = cf.PLASMA_WAVELENGTH_STEPS * delta_x
        f_plasma_hz = C_LIGHT / lambda_plasma
        f_plasma_escalada = f_plasma_hz / DIVISOR_X

        # ¡NUEVO! Cálculo automático de la frecuencia pico
        lambda_peak = cf.STEPS_WAVELENGTH * delta_x
        f_peak_hz = C_LIGHT / lambda_peak
        f_peak_escalada = f_peak_hz / DIVISOR_X

        # Procesamiento del DFT (solo magnitud)
        magnitud = np.abs(datos_dft)

        # ---------------------------------------------------------
        # GENERACIÓN DEL GRÁFICO
        # ---------------------------------------------------------
        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)

        # Trazado de la señal original
        ax.plot(
            frecuencias_escaladas, 
            magnitud, 
            color="black", 
            linewidth=1.5,
            zorder=3
        )

        # Condicional para la línea de plasma
        if MOSTRAR_LINEA_PLASMA:
            ax.axvline(
                x=f_plasma_escalada, 
                color="firebrick",    
                linestyle="--",       
                linewidth=1.3, 
                label=r"Frecuencia de plasma $f_p$",
                zorder=2
            )

        # ¡NUEVO! Condicional para la línea de frecuencia pico
        if MOSTRAR_LINEA_PEAK:
            ax.axvline(
                x=f_peak_escalada, 
                color="royalblue",    # Azul elegante para contrastar con el rojo
                linestyle="--",       
                linewidth=1.3, 
                label=r"Frecuencia pico $f_k$",
                zorder=2
            )

        # Si cualquiera de las dos líneas está activa, muestra la leyenda
        if MOSTRAR_LINEA_PLASMA or MOSTRAR_LINEA_PEAK:
            ax.legend(loc="upper right", frameon=True, edgecolor="none")

        # Configuración dinámica de los ejes
        ax.set_xlabel(f"Frecuencia ({UNIDAD_X})", fontsize=13)
        ax.set_ylabel("Espectro normalizado", fontsize=13)
        
        ax.set_xlim(0, cf.MAX_FREQ / DIVISOR_X)
        ax.set_ylim(0, 1)
        ax.ticklabel_format(style='plain', axis='x', useOffset=False)
        ax.grid(False) 
        plt.tight_layout()
        
        # ---------------------------------------------------------
        # GUARDAR LA IMAGEN
        # ---------------------------------------------------------
        plt.savefig(f"./fotogramas/{NOMBRE_GRAFICO_SALIDA}", dpi=600, bbox_inches='tight')
        print(f"[INFO] Imagen guardada con éxito como: {NOMBRE_GRAFICO_SALIDA}")

        plt.show()

        # Resumen en consola
        print("\n--- Parámetros Físicos Calculados ---")
        print(f"Malla espacial (delta_x) : {delta_x:.4e} m")
        print(f"Paso de tiempo (delta_t) : {delta_t:.4e} s")
        print(f"Tiempo de simulación     : {t_total:.4e} s")
        print(f"Resolución (delta_f)     : {delta_f:.2f} Hz")
        print(f"Frecuencia máxima (lim)  : {cf.MAX_FREQ:.2e} Hz")
        print(f"Frecuencia de Plasma (fp): {f_plasma_hz:.2e} Hz ({f_plasma_escalada:.2f} {UNIDAD_X})")
        print(f"Frecuencia Pico (f0)     : {f_peak_hz:.2e} Hz ({f_peak_escalada:.2f} {UNIDAD_X})") # ¡NUEVO!
        print("-------------------------------------")

    except KeyError:
        print(f"\n[ERROR]: No se encontró el dataset '{DATASET_DFT}' en el archivo HDF5.")
    except AttributeError as e:
        print(f"\n[ERROR]: Falta una variable en tu archivo config.py. Detalles: {e}")
    except FileNotFoundError:
        print(f"\n[ERROR]: El archivo '{ruta_archivo}' no existe en este directorio.")

if __name__ == "__main__":
    graficar_dft_fisica(ARCHIVO_HDF5)