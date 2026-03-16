import sympy as sp
import pyvista as pv

# Define the Golden Ratio
phi = (1 + sp.sqrt(5)) / 2
# Radius of the (0, 1, phi) vertices
rr_sq = 1 + phi ** 2
rr = sp.sqrt(rr_sq)
# Generate the 12 raw vertices (0, +/-1, +/-phi) and its cyclic permutations
raw_vertices = sp.Matrix([
    [0, 1, -rr],
    [0, 1, rr],
    [0, -1, rr],
    [-1, rr, 0],
    [-1, -rr, 0],
    [1, rr, 0],
    [1, -rr, 0],
    [0, -1, -rr],
    [rr, 0, 1],
    [-rr, 0, 1],
    [-rr, 0, -1],
    [rr, 0, -1],
])

icos_points_np = sp.matrix2numpy(raw_vertices, dtype=float)

icos = pv.Icosahedron()
icos.save("loosey.vtp")
icos.points = icos_points_np
icos.save("exacter.vtp")

aa = raw_vertices[0, :].dot(raw_vertices[1, :])
print(sp.acos(aa / rr_sq).simplify())
