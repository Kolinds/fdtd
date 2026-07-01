import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import h5py
import config as cf

# ==========================================
# ESTILO DE PUBLICACIÓN CIENTÍFICA
# ==========================================
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm' 
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

# ==========================================
# CONFIGURACIÓN DE LA CAPTURA Y ENTORNO
# ==========================================
DSET_NAME = cf.EDSET_NAME
FRAME_A_CAPTURAR = 3000
NOMBRE_ARCHIVO_SALIDA = f"captura_onda_t{FRAME_A_CAPTURAR}.pdf"

# --- NUEVAS OPCIONES PARA LA INTERFAZ PLASMA-VACÍO ---
MOSTRAR_INTERFAZ = True  # Cambia a False si no quieres mostrar la separación
NODO_INTERFAZ = 3000      # Cambia este valor al nodo exacto de la interfaz
SOMBREAR_PLASMA = True   # Cambia a False si solo quieres la línea sin el fondo coloreado
# -----------------------------------------------------

# Cargar solo el frame que necesitamos para ahorrar memoria
try:
    with h5py.File(cf.FILE_NAME, "r") as f:
        # Extraemos directamente la fila de datos correspondiente a ese instante
        datos_instante = f[DSET_NAME][FRAME_A_CAPTURAR, :]
        space_size = len(datos_instante)
except Exception as e:
    print(f"[ERROR] No se pudo cargar el archivo HDF5: {e}")
    exit()

# ==========================================
# GENERACIÓN DEL GRÁFICO ESTÁTICO
# ==========================================
fig, ax1 = plt.subplots(figsize=(7, 4.5))

ax1.set_xlim(2500, 3500)
ax1.set_ylim(-1.1, 1.1)
ax1.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])

ax1.set_xlabel(r"Nodo de la malla", fontsize=13)
ax1.set_ylabel(r"Campo eléctrico norm. ($E / E_{max}$)", fontsize=13)

ax1.axhline(0, color='black', linewidth=0.8, linestyle='-.', alpha=0.5)
ax1.grid(True, linestyle=':', alpha=0.6)

# --- DIBUJAR LA INTERFAZ PLASMA-VACÍO ---
if MOSTRAR_INTERFAZ:
    # Línea vertical punteada
    ax1.axvline(x=NODO_INTERFAZ, color='teal', linewidth=1.5, linestyle='--', label=r'Interfaz vacío-plasma')
    
    # Sombreado opcional para identificar visualmente la región del plasma
    if SOMBREAR_PLASMA:
        ax1.axvspan(NODO_INTERFAZ, space_size, color='teal', alpha=0.08)
# ----------------------------------------

# Eje X
space_array = np.arange(space_size)

# Trazado del frame
ax1.plot(space_array, datos_instante, color="firebrick", linewidth=1.5, label=r"$E_z$")

# Añadimos el texto indicando qué frame es
ax1.text(0.03, 0.92, rf'Paso temporal ($n$): {FRAME_A_CAPTURAR}', transform=ax1.transAxes, 
         fontsize=12, bbox=dict(facecolor='white', edgecolor='black', alpha=0.9))

ax1.legend(loc="upper right", framealpha=0.9, edgecolor="black")

# ==========================================
# GUARDAR LA IMAGEN EN ALTA CALIDAD
# ==========================================
plt.savefig(f"./fotogramas/{NOMBRE_ARCHIVO_SALIDA}", dpi=600, bbox_inches='tight')
print(f"Imagen guardada exitosamente como: {NOMBRE_ARCHIVO_SALIDA}")

# Mostrarla en pantalla también para confirmar que quedó bien
plt.show()