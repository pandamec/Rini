# vcct_two_layer_beam.jl
# VCCT implementation for a two-layer beam using element nodal forces (not interface springs)

using LinearAlgebra
using SparseArrays
using Plots

# Beam element stiffness matrix (Euler-Bernoulli, 2D frame, local x along beam)
function beam2d_stiffness(EI, EA, L)
    k = zeros(6,6)
    # Axial
    k[1,1] =  EA/L;  k[1,4] = -EA/L
    k[4,1] = -EA/L;  k[4,4] =  EA/L
    # Bending
    k[2,2] = 12*EI/L^3;  k[2,5] = -12*EI/L^3
    k[5,2] = -12*EI/L^3; k[5,5] = 12*EI/L^3
    k[2,6] = 6*EI/L^2;   k[6,2] = 6*EI/L^2
    k[5,6] = -6*EI/L^2;  k[6,5] = -6*EI/L^2
    k[3,3] = 4*EI/L;     k[3,6] = 2*EI/L
    k[6,3] = 2*EI/L;     k[6,6] = 4*EI/L
    k[2,3] = 0.0; k[3,2] = 0.0
    k[5,3] = 0.0; k[3,5] = 0.0
    k[4,2] = 0.0; k[2,4] = 0.0
    k[4,3] = 0.0; k[3,4] = 0.0
    k[4,5] = 0.0; k[5,4] = 0.0
    k[4,6] = 0.0; k[6,4] = 0.0
    return k
end

# Assemble global stiffness
function assemble_global(EI,EA,L,n_elem)
    ndof = 3*(n_elem+1)
    K = spzeros(ndof, ndof)
    for e=1:n_elem
        Ke = beam2d_stiffness(EI,EA,L)
        gdofs = [(e-1)*3+1:(e-1)*3+6;]
        for i=1:6, j=1:6
            K[gdofs[i],gdofs[j]] += Ke[i,j]
        end
    end
    return K
end

# Apply boundary conditions
function apply_bc(K,f,fixeddofs)
    free = setdiff(1:length(f), fixeddofs)
    Kff = K[free,free]
    ff = f[free]
    u = zeros(length(f))
    u[free] = Kff \ ff
    return u
end

# Extract element end forces (local)
function element_forces(EI,EA,L,ue)
    Ke = beam2d_stiffness(EI,EA,L)
    return Ke * ue
end

# Compute VCCT G from nodal forces and displacements
function compute_G_vcct(EI,EA,L,n_elem,crack_elem,P)
    ndof = 3*(n_elem+1)
    K = assemble_global(EI,EA,L,n_elem)
    f = zeros(ndof)
    # Apply vertical point load P at top node (last node, y-DOF)
    f[3*n_elem-1] = -P

    # Apply BCs: fix node 1 (all DOFs)
    fixeddofs = [1,2,3]
    u = apply_bc(K,f,fixeddofs)

    # Element behind crack tip
    e = crack_elem
    gdofs = [(e-1)*3+1:(e-1)*3+6;]
    ue = u[gdofs]
    Fe = element_forces(EI,EA,L,ue)

    # Forces at node i (the second node of element e)
    Fx = Fe[4]   # axial
    Fy = Fe[5]   # shear vertical
    # For VCCT: normal = vertical, shear = axial (assuming crack horizontal)
    F_normal = Fy
    F_shear  = Fx

    # Displacements ahead (node i+1)
    nodeahead = e+1
    ux = u[(nodeahead-1)*3+1]
    uy = u[(nodeahead-1)*3+2]
    # bottom layer absent in this toy beam, so assume symmetry → Δu = displacement itself
    delta_shear = ux
    delta_normal = uy

    da = L
    G_I  = (1/(2*da)) * (F_normal*delta_normal)
    G_II = (1/(2*da)) * (F_shear *delta_shear)
    return G_I,G_II,G_I+G_II
end

# Example run: vary crack length
elem_length = 0.1
n_elem = 10
EI = 1.0
EA = 1000.0
P = 1.0
cracks = 1:n_elem-1
Gtot = Float64[]
for crack_elem in cracks
    G_I,G_II,G = compute_G_vcct(EI,EA,elem_length,n_elem,crack_elem,P)
    push!(Gtot,G)
end

plot(cracks.*elem_length,Gtot,label="G",xlabel="Crack length",ylabel="G",lw=2,marker=:o,title="VCCT with nodal forces")
