import warp as wp
import warp.fem as fem
import numpy as np
git
wp.init()


# 1. Define the custom Integrand
@fem.integrand
def optimal_transport_loss(
        domain: fem.Domain,
        rho: fem.Field,
        vel: fem.Field,
        alpha: float,
        beta: float
):
    # Evaluate velocity (piecewise constant in cell)
    v = vel(domain)

    # Evaluate density gradient (constant in cell for P1)
    # grad_rho is a 3D vector: [drho/dx, drho/dy, drho/dz]
    grad_rho = fem.grad(rho, domain)

    # Extract components
    drho_dx = grad_rho[0]
    drho_dy = grad_rho[1]
    drho_dz = grad_rho[2]

    # Compute squared velocity norm
    v_norm_sq = v[0] * v[0] + v[1] * v[1]

    # Compute continuity residual (dz/dt + v \cdot \nabla_{xy} rho)
    residual = drho_dz + v[0] * drho_dx + v[1] * drho_dy

    # Return weighted L2 norm squared
    return alpha * v_norm_sq + beta * (residual * residual)


# 2. Setup Mesh and Spaces
# Assume `vertices` (Nx3) and `indices` (Mx4) are your tetrahedral mesh arrays
geo = fem.Tetmesh(vertices=vertices, indices=indices)

# rho is a scalar field on vertices (P1)
space_rho = fem.make_polynomial_space(geo, degree=1)

# vel is a 2D vector field on cells (P0)
space_vel = fem.make_polynomial_space(geo, degree=0, dtype=wp.vec2)

# Create Warp arrays for the unknown degrees of freedom (DOFs)
rho_dof = wp.zeros(space_rho.node_count(), dtype=float, requires_grad=True)
vel_dof = wp.zeros(space_vel.node_count(), dtype=wp.vec2, requires_grad=True)

# Define fields from the DOFs
rho_field = space_rho.make_field(rho_dof)
vel_field = space_vel.make_field(vel_dof)

# 3. Optimization Loop
alpha_weight = 1.0
beta_weight = 1.0
loss = wp.zeros(1, dtype=float, requires_grad=True)

# Create a mask for vertices on z=0 and z=1 to enforce Dirichlet boundaries
# (You would pre-compute this boolean array based on your vertex coordinates)
boundary_mask = wp.array(is_boundary_vertex_array, dtype=bool)


@wp.kernel
def mask_gradients_kernel(grad: wp.array(dtype=float), mask: wp.array(dtype=bool)):
    i = wp.tid()
    if mask[i]:
        grad[i] = 0.0


optimizer_steps = 1000
learning_rate = 0.01

for step in range(optimizer_steps):
    tape = wp.Tape()
    with tape:
        loss.zero_()
        # Integrate the loss over the tetrahedral mesh
        fem.integrate(
            optimal_transport_loss,
            geo,
            fields={"rho": rho_field, "vel": vel_field},
            values={"alpha": alpha_weight, "beta": beta_weight},
            output=loss
        )

    # Compute gradients w.r.t rho_dof and vel_dof
    tape.backward(loss)

    # Zero out gradients on the boundary vertices to lock their values
    wp.launch(
        mask_gradients_kernel,
        dim=space_rho.node_count(),
        inputs=[rho_dof.grad, boundary_mask]
    )

    # Take an optimization step (simple Gradient Descent shown here)
    # For production, export these arrays to an optimizer via DLPack
    rho_grad = rho_dof.grad.numpy()
    vel_grad = vel_dof.grad.numpy()

    # Update DOFs (pseudo-code, you'd update the Warp arrays directly)
    rho_dof_np = rho_dof.numpy() - learning_rate * rho_grad
    vel_dof_np = vel_dof.numpy() - learning_rate * vel_grad

    rho_dof.assign(rho_dof_np)
    vel_dof.assign(vel_dof_np)

    # Clear gradients for the next step
    tape.zero()