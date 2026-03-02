import matplotlib.pyplot as plt
import numpy as np
import h5py

import matplotlib.animation as animation

plt.style.use('dark_background')
fig, ax1 = plt.subplots(figsize=(5, 5), facecolor="#232230",
                 layout='constrained')

ax1.spines['bottom'].set_position(('data', 0))
ax1.set(xlim = (0, 200),ylim = (-1, 1))
ax1.spines["bottom"].set_color("black")
ax1.spines["bottom"].set_linewidth(2.0)
for borde in ("left", "right"):
    ax1.spines[borde].set(color="black", linewidth=2.0)


H5FILE_NAME="wave_data.hdf5"
DATASET_NAME = "elec_fdata"
TOTAL_TIME=500
SPACE_SIZE = 200


space_array = np.linspace(0, SPACE_SIZE, SPACE_SIZE)
wave_line, = ax1.plot(space_array, np.zeros_like(space_array), color = "red")
time_text = ax1.text(0.05, 0.95, 'CUCARACHA', transform=ax1.transAxes, fontsize=12)

with h5py.File(H5FILE_NAME, "r") as f:
        all_data = f[DATASET_NAME][:TOTAL_TIME, :SPACE_SIZE]

def anime_func(frame):
    present_data = all_data[frame, :]
    wave_line.set_ydata(present_data)
    time_text.set_text(f'Frame: {frame}')
    return wave_line, time_text

ani = animation.FuncAnimation(fig=fig, func=anime_func, frames=TOTAL_TIME, interval=3)

plt.show()