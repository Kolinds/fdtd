import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import config as cf

# ==========================================
# CONFIGURACIÓN
# ==========================================
ARCHIVO_HDF5 = "reflec_trans1.hdf5"
DATASET_TRANSMISSION = "transmission"

NOMBRE_GRAFICO_SALIDA = "transmitancia.pdf"

# Escala del eje X
DIVISOR_X = 1e9      # 1=Hz, 1e3=kHz, 1e6=MHz, 1e9=GHz
UNIDAD_X = "GHz"

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


def graficar_transmitancia(ruta_archivo):
    """
    Lee el coeficiente de transmisión complejo t(f)
    almacenado en el HDF5 y grafica la transmitancia:

        T(f) = |t(f)|²
    """

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

        # -------------------------------------------------
        # EJE DE FRECUENCIAS
        # -------------------------------------------------
        frecuencias_hz = np.arange(len(t_complex)) * delta_f
        frecuencias = frecuencias_hz / DIVISOR_X

        # -------------------------------------------------
        # TRANSMITANCIA
        # -------------------------------------------------
        transmitancia = np.abs(t_complex)

        # -------------------------------------------------
        # GRÁFICA
        # -------------------------------------------------
        fig, ax = plt.subplots(figsize=(7, 4.5), dpi=150)

        ax.plot(
            frecuencias,
            transmitancia,
            color="black",
            linewidth=1.5,
            zorder=3
        )

        ax.set_xlabel(f"Frecuencia ({UNIDAD_X})", fontsize=13)
        ax.set_ylabel("Transmitancia", fontsize=13)

        ax.set_xlim((25*delta_f)/DIVISOR_X, (len(t_complex) * delta_f)/DIVISOR_X)
        # ax.set_xlim(0, cf.MAX_FREQ / DIVISOR_X)

        # Ajuste automático del eje Y
        ymax = np.max(transmitancia)

        if ymax <= 1.1:
            ax.set_ylim(0, 1.05)
        else:
            ax.set_ylim(0, 1.05 * ymax)

        ax.ticklabel_format(
            style='plain',
            axis='x',
            useOffset=False
        )

        ax.grid(False)

        plt.tight_layout()

        # -------------------------------------------------
        # GUARDAR FIGURA
        # -------------------------------------------------
        plt.savefig(
            f"./fotogramas/{NOMBRE_GRAFICO_SALIDA}",
            dpi=600,
            bbox_inches='tight'
        )

        print(f"[INFO] Imagen guardada: {NOMBRE_GRAFICO_SALIDA}")

        plt.show()

        # -------------------------------------------------
        # INFORMACIÓN DE CONTROL
        # -------------------------------------------------
        print("\n--- Parámetros FDTD ---")
        print(f"delta_x           = {delta_x:.4e} m")
        print(f"delta_t           = {delta_t:.4e} s")
        print(f"Tiempo total      = {t_total:.4e} s")
        print(f"delta_f           = {delta_f:.4e} Hz")
        print(f"Frecuencia máxima = {cf.MAX_FREQ:.4e} Hz")
        print(f"Número de puntos  = {len(t_complex)}")
        print("-----------------------")

    except KeyError:
        print(
            f"\n[ERROR] No se encontró el dataset "
            f"'{DATASET_TRANSMISSION}' en el archivo HDF5."
        )

    except FileNotFoundError:
        print(
            f"\n[ERROR] El archivo '{ruta_archivo}' no existe."
        )

    except AttributeError as e:
        print(
            f"\n[ERROR] Falta alguna variable en config.py:\n{e}"
        )


if __name__ == "__main__":
    graficar_transmitancia(ARCHIVO_HDF5)