function beam_element(E, I, L)
    k = E * I / L^3 * [
        12   6*L   -12   6*L
        6*L  4*L^2 -6*L  2*L^2
        -12  -6*L   12  -6*L
        6*L  2*L^2 -6*L  4*L^2
    ]
    return k
end

function cohesive_element(K0, L)
    k = K0 * L / 2 * [
        1  0  -1  0
        0  0   0  0
       -1  0   1  0
        0  0   0  0
    ]
    return k
end

function assemble_system(TestSetup,force::Float64,Si,Parylene,Steel,CZM)
    n_dof_steel = length(TestSetup.nodes) * 2
    n_dof_si = (TestSetup.n_elem_layered + 1) * 2
    n_dof = n_dof_steel + n_dof_si
    K = spzeros(n_dof, n_dof)
    F = zeros(n_dof)
    
    I_si = Si.thickness^3 / 12
    I_par = Parylene.thickness^3 / 12
    I_steel = Steel.thickness^3 / 12
    
    for i in 1:TestSetup.n_elem
        idx_steel = [(i-1)*2+1:(i-1)*2+2; (i)*2+1:(i)*2+2] # DoF
        k_steel = beam_element(Steel.E, I_steel, TestSetup.dx)
        K[idx_steel, idx_steel] += k_steel
        
        if i >= TestSetup.start_elem && i <= TestSetup.end_elem
            k_par = beam_element(Parylene.E, I_par, TestSetup.dx)
            K[idx_steel, idx_steel] += k_par
        end
    end
    
    for i in TestSetup.start_elem:TestSetup.end_elem
        idx_steel = [(i-1)*2+1:(i-1)*2+2; (i)*2+1:(i)*2+2]
        idx_si = [n_dof_steel + (i-TestSetup.start_elem)*2+1:n_dof_steel + (i-TestSetup.start_elem)*2+2;
                  n_dof_steel + (i-TestSetup.start_elem+1)*2+1:n_dof_steel + (i-TestSetup.start_elem+1)*2+2]
        
        k_si = beam_element(Si.E, I_si, TestSetup.dx)
        K[idx_si, idx_si] += k_si
        
        k_coh = cohesive_element(CZM.K0, TestSetup.dx)
        K[idx_si, idx_si] += k_coh
        K[idx_steel, idx_steel] += k_coh
        K[idx_si, idx_steel] -= k_coh
        K[idx_steel, idx_si] -= k_coh
    end
    
    mid = div(TestSetup.n_elem, 2) + 1
    F[mid*2-1] = -force
    
    K[1,1] += 1e8  # Left steel support
    K[n_dof_steel-1,n_dof_steel-1] += 1e8  # Right steel support
    # Fix Si rotational DOFs at boundaries of layered region to eliminate rank deficiency
    K[n_dof_steel+2,n_dof_steel+2] += 1e8  # Left Si rotation (29 mm)
    K[n_dof-1,n_dof-1] += 1e8  # Right Si rotation (31 mm)
    
    return K, F
end


function fatigue_degradation!(TestSetup,CZM,K, u, cycles, damage, n_dof_steel)
   damage_hist=zeros(TestSetup.n_elem_layered)
    for i in 1:TestSetup.n_elem_layered
        elem = TestSetup.start_elem + i - 1
        idx_steel = (elem-1)*2 + 1
        idx_si = n_dof_steel + (i-1)*2 + 1
        δ = abs(u[idx_si] - u[idx_steel])
        
        if δ > 0
            damage[i] += CZM.m * (δ / CZM.δ_c) * cycles
            damage[i] = min(damage[i], 1.0)
            push!(damage_hist,CZM.m * (δ / CZM.δ_c) * cycles)

            K_factor = (1 - damage[i])
            idx_steel_full = [(elem-1)*2+1:(elem-1)*2+2; (elem)*2+1:(elem)*2+2]
            idx_si_full = [n_dof_steel + (i-1)*2+1:n_dof_steel + (i-1)*2+2;
                          n_dof_steel + i*2+1:n_dof_steel + i*2+2]
            
            k_coh_old = cohesive_element(CZM.K0, TestSetup.dx)
            k_coh_new = cohesive_element(CZM.K0 * K_factor, TestSetup.dx)
            K[idx_si_full, idx_si_full] += (k_coh_new - k_coh_old)
            K[idx_steel_full, idx_steel_full] += (k_coh_new - k_coh_old)
            K[idx_si_full, idx_steel_full] -= (k_coh_new - k_coh_old)
            K[idx_steel_full, idx_si_full] -= (k_coh_new - k_coh_old)
        end
    end
    return K,damage_hist
end

function simulate_fatigue(TestSetup,max_force::Float64, n_cycles::Int,Si,Parylene,Steel,CZM)
    damage = zeros(TestSetup.n_elem_layered)
    u_history = Vector{Vector{Float64}}()
    n_dof_steel = length(TestSetup.nodes) * 2  # Define locally
    
    K_base, F = assemble_system(TestSetup,max_force,Si,Parylene,Steel,CZM)
    println("Rank of K: ", rank(K_base), " Size: ", size(K_base, 1))
    K = copy(K_base)
    u_initial = K \ F
    println("Initial Max Steel Deflection (μm): ", maximum(abs.(u_initial[1:2:n_dof_steel])) * 1e6)
    println("Initial Max Si Deflection (μm): ", maximum(abs.(u_initial[n_dof_steel+1:2:end])) * 1e6)
    
    damage_history=Vector{Vector{Float64}}()

    for cycle in 1:n_cycles
        u = K \ F
        push!(u_history, copy(u))

        K = copy(K_base)
        K,damage_hist = fatigue_degradation!(TestSetup,CZM,K, u, 1.0, damage, n_dof_steel)
        push!(damage_history,damage_hist)
    end
    
    return u_history, damage, damage_history
end

