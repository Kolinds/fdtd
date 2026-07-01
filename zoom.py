import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import config as cf

# ==========================================
# CONFIGURACIÓN
# ==========================================
ARCHIVOS_HDF5 = ["reflec_trans1.hdf5", "reflec_trans2.hdf5", "reflec_trans3.hdf5"]
ETIQUETAS = ["ADE", "PLRC", "Z-Transform"]
COLORES = ["black", "firebrick", "royalblue", "forestgreen", "darkorchid", "darkslategray"]

DATASET_TRANSMISSION = "reflection"

# --- INTERRUPTOR DE ESCALA ---
ESCALA_DB = False  # True para escala logarítmica (dB), False para escala lineal
# -----------------------------

# --- CONFIGURACIÓN DEL ZOOM (INSET) ---
SELECCION_INTERACTIVA = False  # True: Elige la región haciendo 2 clics. False: Usa los límites de abajo.

# Si SELECCION_INTERACTIVA = False, se usarán estos límites fijos para el zoom:
ZOOM_X = (0.209, 0.213)          
ZOOM_Y = (0.35, 0.37)            

# Posición del gráfico flotante (Inset) dentro de la gráfica principal.
# Formato: [x_origen, y_origen, ancho, alto] en fracciones del 0 al 1.
# Ejemplo: [0.55, 0.55, 0.4, 0.4] lo coloca en la esquina superior derecha.
UBICACION_INSET = [0.55, 0.50, 0.42, 0.42]
# -------------------------------

# --- MARCADOR DE FRECUENCIA DE PLASMA ---
MOSTRAR_FREQ_PLASMA = True  
# ----------------------------------------

NOMBRE_GRAFICO_BASE = "transmitancia_multiples_casos_inset"
DIVISOR_X = 1e6      
UNIDAD_X = "MHz"
C_LIGHT = 299792458.0
# ==========================================

# ==========================================
# ESTILO DE PUBLICACIÓN
# ==========================================
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.size'] = 11
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
    Grafica la transmisión completa y añade un gráfico Inset (flotante) 
    conducido por líneas de expansión.
    """
    # Volvemos a un solo lienzo principal
    fig, ax1 = plt.subplots(figsize=(8, 5), dpi=150)
    
    # Creamos el eje flotante (Inset) dentro de ax1
    ax_inset = ax1.inset_axes(UBICACION_INSET)
    
    longitud_puntos = 0
    delta_f_control = 0
    delta_x_control = 0
    delta_t_control = 0
    t_total_control = 0

    # Iteramos sobre cada archivo para rellenar las curvas
    for i, ruta_archivo in enumerate(lista_archivos):
        try:
            with h5py.File(ruta_archivo, "r") as f:
                t_complex = f[DATASET_TRANSMISSION][:]

            lambda_min = C_LIGHT / cf.MAX_FREQ
            delta_x = lambda_min / 20.0
            delta_t = cf.COURANT * (delta_x / C_LIGHT)
            t_total = cf.TOTAL_TIME * delta_t
            delta_f = 1.0 / t_total

            longitud_puntos = len(t_complex)
            delta_f_control = delta_f
            delta_x_control = delta_x
            delta_t_control = delta_t
            t_total_control = t_total

            frecuencias_hz = np.arange(len(t_complex)) * delta_f
            frecuencias = frecuencias_hz / DIVISOR_X
            
            magnitud_lineal = np.abs(t_complex)
            
            if ESCALA_DB:
                datos_y = 20 * np.log10(magnitud_lineal + 1e-12)
            else:
                datos_y = magnitud_lineal

            # Trazamos en la gráfica principal y en la flotante
            ax1.plot(frecuencias, datos_y, color=COLORES[i], label=lista_etiquetas[i], linewidth=1.5, zorder=3)
            ax_inset.plot(frecuencias, datos_y, color=COLORES[i], linewidth=1.5, zorder=3)

        except (KeyError, FileNotFoundError, AttributeError) as e:
            print(f"[ERROR] Error al procesar el archivo {ruta_archivo}: {e}")
            return

    # --- CÁLCULO DINÁMICO DE FRECUENCIA DE PLASMA ---
    if MOSTRAR_FREQ_PLASMA:
        try:
            lambda_plasma = cf.PLASMA_WAVELENGTH_STEPS * delta_x_control
            freq_plasma_hz = C_LIGHT / lambda_plasma
            freq_plasma_grafico = freq_plasma_hz / DIVISOR_X

            ax1.axvline(x=freq_plasma_grafico, color='darkorange', linewidth=1.5, linestyle='--', label=r'$f_p$', zorder=4)
            ax_inset.axvline(x=freq_plasma_grafico, color='darkorange', linewidth=1.5, linestyle='--', zorder=4)
        except AttributeError:
            print("[ERROR] No se encontró 'PLASMA_WAVELENGTH_STEPS' en config.py")

    # --- CONFIGURACIÓN GRÁFICA PRINCIPAL (ax1) ---
    ax1.set_xlabel(f"Frecuencia ({UNIDAD_X})", fontsize=12)
    ax1.grid(True, linestyle=":", alpha=0.5, color="gray")
    ax1.ticklabel_format(style='plain', axis='x', useOffset=False)

    if ESCALA_DB:
        ax1.set_ylabel("Coef. transmisión (dB)", fontsize=12)
        ax1.set_ylim(-80, 2.0)
    else:
        ax1.set_ylabel("Coef. reflexión (u.a.)", fontsize=12)
        ax1.set_ylim(-0.05, 1.1)
        
    ax1.set_xlim(0.1, 0.7)
    ax1.legend(loc="lower left", frameon=True, edgecolor="none", fontsize=10)

    # --- GESTIÓN DEL ZOOM (INTERACTIVO O MANUAL) ---
    limites_x, limites_y = ZOOM_X, ZOOM_Y

    if SELECCION_INTERACTIVA:
        print("\n[INTERACTIVO] Haz DOS CLICS en el gráfico para definir el área de ZOOM.")
        puntos = plt.ginput(2, timeout=0) 
        
        if len(puntos) == 2:
            x_coords = [puntos[0][0], puntos[1][0]]
            y_coords = [puntos[0][1], puntos[1][1]]
            limites_x = (min(x_coords), max(x_coords))
            limites_y = (min(y_coords), max(y_coords))
        else:
            print("[ADVERTENCIA] Selección incompleta. Usando valores por defecto.")

    # --- CONFIGURACIÓN GRÁFICO FLOTANTE (ax_inset) ---
    ax_inset.set_xlim(limites_x)
    ax_inset.set_ylim(limites_y)
    ax_inset.grid(True, linestyle=":", alpha=0.5, color="gray")
    ax_inset.tick_params(axis='both', which='major', labelsize=9) # Números más pequeños para que quepan bien
    
    # Fondo blanco opaco para que las curvas principales no interfieran visualmente
    ax_inset.set_facecolor('white')

    # =========================================================================
    # EFECTO DE ZOOM: CUADRO CON LÍNEAS CONECTORAS
    # =========================================================================
    # indicate_inset_zoom ahora conecta ax1 con ax_inset
    rectangulo_indicador, lineas_conectoras = ax1.indicate_inset_zoom(
        ax_inset, 
        edgecolor="dimgray",    
        linewidth=1.5,          
        alpha=1.0,              
        linestyle="-",          
        zorder=5                
    )
    
    for conector in lineas_conectoras:
        conector.set_visible(True)              
        conector.set_color("firebrick")         
        conector.set_linewidth(1.8)             
        conector.set_alpha(0.8)                 
        conector.set_linestyle((0, (4, 3)))     
    # =========================================================================

    plt.tight_layout()

    # Guardar y mostrar
    nombre_salida = f"{NOMBRE_GRAFICO_BASE}_db.pdf" if ESCALA_DB else f"{NOMBRE_GRAFICO_BASE}_lineal.pdf"
    plt.savefig(f"./fotogramas/{nombre_salida}", dpi=600, bbox_inches='tight')
    print(f"[INFO] Imagen Inset guardada exitosamente como: ./fotogramas/{nombre_salida}")
    plt.show()


if __name__ == "__main__":
    graficar_transmitancias_multiples(ARCHIVOS_HDF5, ETIQUETAS)