import numpy as np


def ezfield(current_time, location, time_delay, loc_delay, width, courant):
    return np.exp((-(current_time + time_delay + loc_delay - location / courant)**2)/width)