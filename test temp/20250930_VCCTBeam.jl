using LinearAlgebra

# Parameters
E = 1e9  # Young's modulus (MPa)
h = 3.0  # Thickness of each arm (mm)
b = 20.0  # Width (mm)
L = 150.0  # Total length (mm)
a0 = 50.0  # Initial crack length (mm)
Gc = 0.3  # Critical energy release rate (N/mm)
delta_max = 5.0  # Maximum total opening displacement (mm)
n_inc = 50  # Number of displacement increments
n_elem = 100  # Number of elements along the length
dx = L / n_elem
I = b * h^3 / 12  # Moment of inertia
EI = E * I
alpha = 1e10  # Penalty parameter for bonding constraints

# Number of nodes
n_nodes = n_elem + 1

# Total degrees of freedom (v_upper, theta_upper, v_lower, theta_lower per node)
n_dof = 4 * n_nodes

# Beam element stiffness matrix
function beam_k(le, EI)
    k = zeros(4, 4)
    k[1, 1] = 12 * EI / le^3
    k[1, 2] = 6 * EI / le^2
    k[1, 3] = -12 * EI / le^3
    k[1, 4] = 6 * EI / le^2
    k[2, 1] = 6 * EI / le^2
    k[2, 2] = 4 * EI / le
    k[2, 3] = -6 * EI / le^2
    k[2, 4] = 2 * EI / le
    k[3, 1] = -12 * EI / le^3
    k[3, 2] = -6 * EI / le^2
    k[3, 3] = 12 * EI / le^3
    k[3, 4] = -6 * EI / le^2
    k[4, 1] = 6 * EI / le^2
    k[4, 2] = 2 * EI / le
    k[4, 3] = -6 * EI / le^2
    k[4, 4] = 4 * EI / le
    return k
end

# Assemble global stiffness matrix
function assemble_K(current_p)
    K = zeros(n_dof, n_dof)
    le = dx
    k = beam_k(le, EI)
    for e in 1:n_elem
        # Upper beam element
        i1 = e
        i2 = e + 1
        d1 = 4 * (i1 - 1) + 1  # v_upper i1
        d2 = d1 + 1  # theta_upper i1
        d3 = 4 * (i2 - 1) + 1  # v_upper i2
        d4 = d3 + 1  # theta_upper i2
        dofs = [d1, d2, d3, d4]
        K[dofs, dofs] .+= k

        # Lower beam element
        d1 = 4 * (i1 - 1) + 3  # v_lower i1
        d2 = d1 + 1  # theta_lower i1
        d3 = 4 * (i2 - 1) + 3  # v_lower i2
        d4 = d3 + 1  # theta_lower i2
        dofs = [d1, d2, d3, d4]
        K[dofs, dofs] .+= k
    end

    # Add penalty constraints for bonded nodes (from current_p+1 to n_nodes)
    for i in (current_p + 1):n_nodes
        duv = 4 * (i - 1) + 1  # v_upper
        dut = duv + 1  # theta_upper
        dlv = duv + 2  # v_lower
        dlt = dut + 2  # theta_lower

        # Penalty for v_upper - v_lower = 0
        K[duv, duv] += alpha
        K[duv, dlv] -= alpha
        K[dlv, duv] -= alpha
        K[dlv, dlv] += alpha

        # Penalty for theta_upper - theta_lower = 0
        K[dut, dut] += alpha
        K[dut, dlt] -= alpha
        K[dlt, dut] -= alpha
        K[dlt, dlt] += alpha
    end
    return K
end

# Apply boundary conditions (modify K and F for Dirichlet BCs)
function apply_bc!(K, F, delta)
    # v_upper at node 1 = delta/2
    d = 1
    K[d, :] .= 0
    K[:, d] .= 0
    K[d, d] = 1
    F[d] = delta / 2

    # v_lower at node 1 = -delta/2
    d = 3
    K[d, :] .= 0
    K[:, d] .= 0
    K[d, d] = 1
    F[d] = -delta / 2

    # v_lower at last node = 0
    d = 4 * (n_nodes - 1) + 3
    K[d, :] .= 0
    K[:, d] .= 0
    K[d, d] = 1
    F[d] = 0

    # theta_lower at last node = 0
    d = 4 * (n_nodes - 1) + 4
    K[d, :] .= 0
    K[:, d] .= 0
    K[d, d] = 1
    F[d] = 0
end

# Compute energy release rate G using VCCT
function compute_G(u, current_p)
    if current_p >= n_nodes - 1
        return 0.0  # No more propagation
    end
    behind_i = current_p
    tip_i = current_p + 1

    # Opening displacement at node behind tip
    duv = 4 * (behind_i - 1) + 1
    dlv = 4 * (behind_i - 1) + 3
    delta_v = u[duv] - u[dlv]

    # Displacement difference at tip (small due to penalty)
    duv_tip = 4 * (tip_i - 1) + 1
    dlv_tip = 4 * (tip_i - 1) + 3
    delta_v_tip = u[duv_tip] - u[dlv_tip]

    # Force at tip connection
    F_y = alpha * delta_v_tip

    # Energy release rate (mode I)
    G = F_y * delta_v / (2 * b * dx)
    return G
end

# Main simulation loop
function simulate_crack_propagation()
    current_p = Int(floor(a0 / dx)) + 1  # Ensure integer
    println("Starting simulation with initial crack length a0 = $a0 mm")

    for inc in 1:n_inc
        delta = (inc / n_inc) * delta_max
        propagated = true
        G = 0.0
        u = zeros(n_dof)
        while propagated
            if current_p >= n_nodes
                println("Crack has propagated beyond the model length.")
                break
            end
            K = assemble_K(current_p)
            F = zeros(n_dof)
            apply_bc!(K, F, delta)
            u = K \ F
            G = compute_G(u, current_p)
            if G > Gc
                current_p += 1
                propagated = (current_p < n_nodes)
                if !propagated
                    println("Crack has propagated to the end.")
                end
            else
                propagated = false
            end
        end

        a = current_p * dx
        # Compute load (reaction at upper loading point)
        K_orig = assemble_K(current_p)
        d_u = 1
        P = dot(K_orig[d_u, :], u)

        println("Increment $inc: delta = $delta mm, a = $a mm, G = $G N/mm, Load P = $P N")
    end
    println("Simulation complete.")
end

# Run the simulation
simulate_crack_propagation()
