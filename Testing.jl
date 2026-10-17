using LinearAlgebra
import PyPlot as plt

m, n = 8, 9

θ = acos(((n^2 + 4*n*m + m^2) / (n^2 + n*m + m^2)) / 2)

# a = 1 angstroms 
a1 = Vector([1, 0])
a2 = Vector([1 / 2, sqrt(3) / 2])

L1 = n * a1 + m * a2
L2 = -m * a1 + (n + m) * a2

latticePoints_1 = [n * a1 + m * a2 for n in -50:50, m in -50:50]
x_1 = [v[1] for v in latticePoints_1]
y_1 = [v[2] for v in latticePoints_1]

R_θ = [cos(θ) -sin(θ);
       sin(θ) cos(θ)]

a1_p = R_θ * a1
a2_p = R_θ * a2

latticePoints_2 = [n * a1_p + m * a2_p for n in -50:50, m in -50:50]
x_2 = [v[1] for v in latticePoints_2]
y_2 = [v[2] for v in latticePoints_2]

plt.scatter(x_1, y_1, marker=".")
plt.scatter(x_2, y_2, marker=".")

L1x = L1[1]
L1y = L1[2]
L2x = L2[1]
L2y = L2[2]

unitCellx = [(v*L1x + u*L2x) / 100 for v in 0:100, u in 0:100]
unitCelly = [(v*L1y + u*L2y) / 100 for v in 0:100, u in 0:100]

plt.scatter(unitCellx, unitCelly, marker=".", alpha=0.1)

plt.show()

# plan: do real space over N1xN2 moire unit cells, do something with dictionaries storing hopping terms, use multiprocessing, plot bands
