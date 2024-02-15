import torch as th
from math import cos, sin, pi
print(th.__version__)
N = 64
a = th.zeros(N, N, dtype=th.cfloat)
for i in range(N):
    for j in range(N):
        a[i, j] = cos((-2 * pi * i * j) / N) + 1.j * sin((-2 * pi * i * j) / N)

b = th.zeros(1, N, dtype=th.cfloat)
for i in range(N):
    b[0, i] = cos(-2 * pi * (i % 3) / 3)  + 1.j * sin(-2 * pi * (i % 3) / 3)

c = b @ a
print(c)