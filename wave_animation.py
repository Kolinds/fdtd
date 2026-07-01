import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib as mpl
import numpy as np
import h5py
import config as cf

# ==========================================
# ESTILO DE PUBLICACIÓN CIENTÍFICA (TFG / Papers)
# ==========================================
mpl.rcParams['font.family'] = 'serif'
# Forzamos la fuente matemática "Computer Modern" (la estándar de LaTeX)
mpl.rcParams['mathtext.fontset'] = 'cm' 
mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.linewidth'] = 1.2
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0

FILE_NAME = "plasma_data.hdf5"
DSET_NAME = cf.EDSET_NAME

# Cargar datos con manejo de errores básico
try:
    with h5py.File(FILE_NAME, "r") as f:
        all_data = f[DSET_NAME][:cf.TOTAL_TIME, :]
        space_size = all_data.shape[1]
except Exception as e:
    print(f"[ERROR] No se pudo cargar el archivo HDF5: {e}")
    exit()

# ==========================================
# CONFIGURACIÓN VISUAL DEL GRÁFICO
# ==========================================
# dpi=150 asegura alta resolución si decides guardar un frame estático
fig, ax1 = plt.subplots(figsize=(7, 4.5), dpi=150, layout='constrained')

ax1.set_xlim(0, space_size)
ax1.set_ylim(-1.1, 1.1)

# Forzamos marcas (ticks) exactas y simétricas para destacar la normalización
ax1.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])

# Etiquetas usando cadenas crudas (r"...") para renderizar LaTeX
ax1.set_xlabel(r"Nodo de la malla", fontsize=13)
ax1.set_ylabel(r"Campo eléctrico normalizado ($E / E_{max}$)", fontsize=13)

# Línea de referencia en el cero (estilo dash-dot, muy clásico en ingeniería)
ax1.axhline(0, color='black', linewidth=0.8, linestyle='-.', alpha=0.5)

# Cuadrícula muy sutil para ayudar a leer la amplitud de los picos
ax1.grid(True, linestyle=':', alpha=0.6)

# ==========================================
# ANIMACIÓN
# ==========================================
# np.arange es ligeramente más eficiente y exacto aquí que linspace
space_array = np.arange(space_size)

# Trazado de la onda usando un rojo oscuro formal ("firebrick") o azul ("navy")
wave_line, = ax1.plot(space_array, np.zeros_like(space_array), 
                      color="firebrick", linewidth=1.5, label=r"$E_z$")

# Caja de texto: Fondo blanco sin bordes bruscos para no distraer
time_text = ax1.text(0.03, 0.92, '', transform=ax1.transAxes, fontsize=12,
                     bbox=dict(facecolor='white', edgecolor='none', alpha=0.85))

# Leyenda (opcional pero le da un toque muy técnico)
ax1.legend(loc="upper right", framealpha=0.9, edgecolor="black")

starting_frame = 0
frame_interval = 100

def anime_func(frame):
    real_frame = (frame * frame_interval) + starting_frame
    
    if real_frame >= cf.TOTAL_TIME:
        real_frame = cf.TOTAL_TIME - 1
        
    present_data = all_data[real_frame, :]
    wave_line.set_ydata(present_data)
    
    # Texto renderizado con la variable temporal n en formato matemático
    time_text.set_text(r'Paso temporal ($n$): ' + str(real_frame))
    
    return wave_line, time_text

frames_totales = (cf.TOTAL_TIME - starting_frame) // frame_interval

ani = animation.FuncAnimation(fig=fig, func=anime_func, frames=frames_totales, 
                              interval=100, blit=True)


plt.show()