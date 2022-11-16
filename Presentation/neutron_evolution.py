import re

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

from entropyOutput import shannonEntropy


def read_generation_positions_from_file(file):
    generations = []

    # https://stackoverflow.com/questions/4703390/how-to-extract-a-floating-number-from-a-string
    floating_point_rx = re.compile(
        "[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?", re.VERBOSE
    )

    with open(file, "r") as file:
        next_gen = []
        for line in file:
            line = line.strip()
            if line == "":
                generations.append(next_gen)
                next_gen = []
                continue

            numbers_on_line = list(map(float, floating_point_rx.findall(line)))
            next_gen.append(numbers_on_line[0:2])

            # Limit animation length
            if len(generations) >= 200:
                break
    return np.array(generations)


positions = read_generation_positions_from_file("positions_large_high_pop.txt")

fig = plt.figure()
zeros = np.zeros(len(positions[0, :, 0]))

(line,) = plt.plot([], color="orange", label="Entropy")
scatter = plt.scatter(zeros, zeros, s=0.1, label="Fission Events")
plt.legend()

xdata, ydata = [], []


def animate_plot(animation_data):
    data, frame = animation_data
    scatter.set_offsets(positions[frame, :, :])

    if frame == 0:
        xdata.clear()
        ydata.clear()

    # TODO: nicer parameter generation, this is rather adhoc
    xdata.append(frame * 4 - 450)
    ydata.append(data * 70 - 550)
    line.set_data(xdata, ydata)


anim = animation.FuncAnimation(
    fig, animate_plot, frames=zip(shannonEntropy, range(len(positions))), save_count=200
)

plt.axis("square")
plt.xlim([-500, 500])
plt.ylim([-500, 500])
anim.save("evolution.mp4")
plt.show()
