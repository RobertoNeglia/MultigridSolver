import numpy as np

N = 8
# A = 49 * 49

# IhH = 4 x 49
# r = 49 x 1
# IhH * r = (4 x 49) * (49 * 1) = (4 x 1)

print("Building interpolator")
# Interpolator (assume N even)
NI, MI = (N - 2) // 2, N - 1 # NI = 3, MI = 7
IhH = np.zeros((NI * NI, MI * MI)) # rows = 9, cols = 49

for il in range(NI):
    i = 1 + il * 2
    for jl in range(NI):
        ijl = il + jl * NI
        j = 1 + jl * 2

        IhH[ijl, i + j * MI] = 1.0 / 4.0 # 0.25 on the middle element
        IhH[ijl, i - 1 + j * MI] = 1.0 / 8.0 # 0.125 on the righ - left - up - down
        IhH[ijl, i + 1 + j * MI] = 1.0 / 8.0
        IhH[ijl, i + (j - 1) * MI] = 1.0 / 8.0
        IhH[ijl, i + (j + 1) * MI] = 1.0 / 8.0
        IhH[ijl, i + 1 + (j - 1) * MI] = 1.0 / 16.0 # 0.0625 on the diagonals
        IhH[ijl, i + 1 + (j + 1) * MI] = 1.0 / 16.0
        IhH[ijl, i - 1 + (j - 1) * MI] = 1.0 / 16.0
        IhH[ijl, i - 1 + (j + 1) * MI] = 1.0 / 16.0

print(IhH)
print(IhH.shape)