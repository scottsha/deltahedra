import pyvista as pv
import numpy as np


def make_rot(point: np.ndarray, axis: np.ndarray) -> np.ndarray:
    shift = np.identity(4)
    shift[0, 3] = point[0]
    shift[1, 3] = point[1]
    shift[2, 3] = point[2]
    rot = np.identity(4)
    rot[0:3, 0:3] += -2 * np.outer(axis, axis)
    back_shift = np.identity(4)
    back_shift[0, 3] = -point[0]
    back_shift[1, 3] = -point[1]
    back_shift[2, 3] = -point[2]
    return shift @ rot @ back_shift


icos = pv.Icosahedron()


def make_pent_flip(link_vert: int) -> np.ndarray:
    r_dir = icos.points[link_vert, :]
    r_dir = r_dir / np.linalg.norm(r_dir)
    flip_link = icos.point_neighbors(link_vert)
    r_center = icos.points[flip_link, :].mean(axis=0)
    rr = make_rot(r_center, r_dir)
    return rr


x = icos.points[:, 0]
y = icos.points[:, 1]
icos.point_data["og"] = np.square(x - 0.2) - y
icos.point_data["id"] = np.array(range(12))


r0 = np.identity(4)
r0[0, 0] = -1.0

phi = (1 + np.sqrt(5)) / 2

phi_aa = (1 - phi) / 2
phi_bb = phi/2

r1 = np.array([
    [phi_aa, -phi_bb, -0.5, 0.0],
    [-phi_bb, 0.5, phi_aa, 0.0],
    [-0.5, phi_aa, phi_bb, 0.0],
    [0.0, 0.0, 0.0, 1.0]
])

r2 = np.identity(4)
r2[1, 1] = -1.0

r3 = make_pent_flip(0)

r5 =

tri_flip_face = 0
r4_dir = icos.face_normals[tri_flip_face, :]
face_0 = icos.regular_faces[tri_flip_face]
print(face_0)
r4_center = icos.points[face_0, :].mean(axis=0)
r4 = make_rot(r4_center, r4_dir)
print(r4)



icos.save("icos.vtp")
icos.transform(r0, inplace=False).save("icos_r0.vtp")
icos.transform(r1, inplace=False).save("icos_r1.vtp")
icos.transform(r2, inplace=False).save("icos_r2.vtp")
icos.transform(r3, inplace=False).save("icos_r3.vtp")
icos.transform(r4, inplace=False).save("icos_r4.vtp")


r0_dir = np.array([1.0, 0.0, 0.0])
r1_dir = np.array([phi_bb, 0.5, -phi_aa])
r2_dir = np.array([0.0, 1.0, 0.0])
dirs = [r0_dir, r1_dir, r2_dir, r3_dir, r4_dir]

print()

coxeter = np.zeros(shape=(len(dirs), len(dirs)))
for ii, dir in enumerate(dirs):
    for jj in range(ii):
        theta = np.arccos(np.dot(dir, dirs[jj]))
        coxeter[jj, ii] = np.pi / theta

print(coxeter)