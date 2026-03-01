import matplotlib.pyplot as plt
import numpy as np

import matplotlib.animation as animation

fig = plt.figure(figsize=(5, 5), facecolor='lightskyblue',
                 layout='constrained')
fig.suptitle('Figure')
ax = fig.add_subplot()
ax.set_title('Axes', loc='left', fontstyle='oblique', fontsize='medium')
plt.show()