import matplotlib.pyplot as plt
import numpy as np
import h5py
import matplotlib.animation as animation
import config as cf


#Visual
plt.style.use('dark_background')
fig, ax1 = plt.subplots(figsize=(5, 5), facecolor="#232230",
                 layout='constrained')

ax1.spines['bottom'].set_position(('data', 0))
ax1.set(xlim = (0, cf.SPACE_SIZE),ylim = (-1, 1))
ax1.spines["bottom"].set_color("black")
ax1.spines["bottom"].set_linewidth(2.0)
for borde in ("left", "right"):
    ax1.spines[borde].set(color="black", linewidth=2.0)


#Animacion
space_array = np.linspace(0, cf.SPACE_SIZE, cf.SPACE_SIZE)
wave_line, = ax1.plot(space_array, np.zeros_like(space_array), color = "red")
time_text = ax1.text(0.05, 0.95, 'CUCARACHA', transform=ax1.transAxes, fontsize=12)

with h5py.File(cf.FILE_NAME, "r") as f:
        all_data = f[cf.EDSET_NAME][:cf.TOTAL_TIME, :cf.SPACE_SIZE]

def anime_func(frame):
    frame = frame * 2
    present_data = all_data[frame, :]
    wave_line.set_ydata(present_data)
    time_text.set_text(f'Frame: {frame}')
    return wave_line, time_text

ani = animation.FuncAnimation(fig=fig, func=anime_func, frames=cf.TOTAL_TIME//2, interval = 20)

plt.show()