import sympy as sp
import pyvista as pv


def make_rot(point: sp.Matrix, axis: sp.Matrix) -> sp.Matrix:
    shift = sp.eye(4)
    shift[0, 3] = point[0]
    shift[1, 3] = point[1]
    shift[2, 3] = point[2]
    rot = sp.eye(4)
    a_norm = axis.dot(axis)
    rot[0:3, 0:3] += -2 * axis.transpose() * axis / a_norm
    back_shift = sp.eye(4)
    back_shift[0, 3] = -point[0]
    back_shift[1, 3] = -point[1]
    back_shift[2, 3] = -point[2]
    stack = shift * rot * back_shift
    return stack


# Define the Golden Ratio
phi = (1 + sp.sqrt(5)) / 2
# Radius of the (0, 1, phi) vertices
rr_sq = 1 + phi ** 2
rr = sp.sqrt(rr_sq)
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


r0 = sp.eye(4)
r0[0, 0] = -1

phi_aa = (1 - phi) / 2
phi_bb = phi/2
r1 = sp.Matrix([
    [phi_aa, -phi_bb, -1/2, 0],
    [-phi_bb, 1/2, phi_aa, 0],
    [-1/2, phi_aa, phi_bb, 0],
    [0, 0, 0, 1]
])

r2 = sp.eye(4)
r2[1, 1] = -1

flip_vert = 3
r3_dir = icos.points[flip_vert, :]
flip_link = icos.point_neighbors(flip_vert)
r3_center = icos.points[flip_link, :].mean(axis=0)
r3 = make_rot(r3_center, r3_dir)
print(r3)


icos.save("icos.vtp")
icos.transform(sp.matrix2numpy(r0, dtype=float), inplace=False).save("icos_r0.vtp")
icos.transform(sp.matrix2numpy(r1, dtype=float), inplace=False).save("icos_r1.vtp")
icos.transform(sp.matrix2numpy(r2, dtype=float), inplace=False).save("icos_r2.vtp")
icos.transform(sp.matrix2numpy(r3, dtype=float), inplace=False).save("icos_r3.vtp")
# icos.transform(sp.matrix2numpy(r4, dtype=float), inplace=False).save("icos_r3.vtp")
# icos.transform(r4, inplace=False).save("icos_r4.vtp")
