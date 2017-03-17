#!/usr/bin/env python3
import matplotlib.pyplot as plt
from math import sqrt
from math import log

dx = [1/sqrt(16), 1/sqrt(64), 1/sqrt(256), 1/sqrt(1024)]
#dx = [1/sqrt(16), 1/sqrt(64), 1/sqrt(256)]

#rl2_euler = [0.00325969, 0.00114543, 0.000456617, 0.000206493]
rl2_euler = [0.0032597, 0.00114543, 0.000456618, 0.000206494]

rl2_ns = [0.00312787, 0.00134238, 0.00150916, 0.00285177]

rul2_euler = [0.0190524, 0.00670444, 0.00289378, 0.00148658]

print(log(rl2_euler[0]/rl2_euler[1])/log(dx[0]/dx[1]))
print(log(rl2_ns[0]/rl2_ns[1])/log(dx[0]/dx[1]))

plt.loglog(dx, rl2_euler, dx, rl2_ns)
plt.grid(True,which="both")
plt.show()

