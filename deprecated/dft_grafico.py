import h5py
import matplotlib.pyplot as plt
import numpy as np
import config as cf  # Importamos tu archivo de configuración

# ==========================================
# CONFIGURACIÓN: Rutas y Datasets
# ==========================================
ARCHIVO_HDF5 = "wave_data.hdf5"
DATASET_DFT = "/Probes/Probe1"  

# Constante de la luz en el vacío (m/s)
C_LIGHT = 299792458.0 
# ==========================================

def graficar_dft_fisica(ruta_archivo):
    """Lee la DFT y grafica la respuesta en frecuencia real usando datos físicos."""
    try:
        with h5py.File(ruta_archivo, "r") as f:
            datos_dft = f[DATASET_DFT][:]

        # ---------------------------------------------------------
        # FÍSICA Y MATEMÁTICAS: Cálculo de resoluciones (1D por defecto)
        # ---------------------------------------------------------
        # 1. Longitud de onda mínima y mallado espacial (10 puntos por onda)
        lambda_min = C_LIGHT / cf.MAX_FREQ
        delta_x = lambda_min / 20.0

        # 2. Límite de Courant y Paso temporal (Delta t)
        # Nota: Si tu simulación es 2D o 3D, debes dividir C_LIGHT por sqrt(2) o sqrt(3)
        delta_t = cf.COURANT * (delta_x / C_LIGHT)

        # 3. Resolución frecuencial (f_min y delta_f)
        # Asumimos que cf.total_time contiene el número total de iteraciones (pasos)
        t_total = cf.TOTAL_TIME * delta_t
        f_min = 1.0 / t_total
        delta_f = f_min

        # 4. Generar el vector de frecuencias reales en Hertz
        # El índice 0 es 0 Hz (DC), el índice 1 es f_min, etc.
        frecuencias_hz = np.arange(len(datos_dft)) * delta_f
        # ---------------------------------------------------------

        # Procesamiento del DFT
        magnitud = np.abs(datos_dft)
        fase = np.angle(datos_dft, deg=True)

        # Configuración de la gráfica
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # --- Gráfica superior: Magnitud ---
        ax1.plot(
            frecuencias_hz, magnitud, color="darkcyan", linewidth=2, label="Amplitud"
        )
        ax1.set_ylabel("Magnitud / Amplitud")
        ax1.set_title(
            f"Espectro Físico de la Señal ($\Delta t$ = {delta_t:.2e} s)",
            fontsize=14,
            fontweight="bold",
        )
        ax1.grid(True, linestyle="--", alpha=0.5)
        ax1.legend()

        # --- Gráfica inferior: Fase ---
        ax2.plot(frecuencias_hz, fase, color="crimson", linewidth=1.5, label="Fase")
        ax2.set_xlabel("Frecuencia (Hz)")
        ax2.set_ylabel("Fase (Grados)")
        ax2.set_ylim(-185, 185)
        ax2.set_yticks([-180, -90, 0, 90, 180])
        ax2.grid(True, linestyle="--", alpha=0.5)
        ax2.legend()
        
        # Limitar el eje X para que solo muestre hasta f_max (o f_min como inicio visible)
        # Empezamos en 0 (componente continua DC) hasta la frecuencia máxima de interés.
        ax2.set_xlim(0, cf.MAX_FREQ)

        plt.tight_layout()
        plt.show()

        # Resumen en consola para verificación rápida
        print("--- Parámetros Físicos Calculados ---")
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
        print("Asegúrate de tener definidas: f_max, courant y total_time en config.py.")
    except FileNotFoundError:
        print(f"\n[ERROR]: El archivo '{ruta_archivo}' no existe en este directorio.")

if __name__ == "__main__":
    graficar_dft_fisica(ARCHIVO_HDF5)