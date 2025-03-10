function beam_element(E, I, L)
    k = E * I / L^3 * [
        12   6*L   -12   6*L
        6*L  4*L^2 -6*L  2*L^2
        -12  -6*L   12  -6*L
        6*L  2*L^2 -6*L  4*L^2
    ]
    return k
end

function cohesive_element(K0, L,w)
    k = K0 * w*L / 2 * [
        1  0  -1  0
        0  0   0  0
       -1  0   1  0
        0  0   0  0
    ]
    return k
end

function assemble_system(TestSetup,force::Float64,Si,Parylene,Steel,CZM)
    n_dof_steel = length(TestSetup.nodes) * 2
    n_dof_si    = (TestSetup.n_elem_layered + 1) * 2
    n_dof       = n_dof_steel + n_dof_si
    K           = spzeros(n_dof, n_dof)
    F           = zeros(n_dof)
    
    I_si        =  Si.thickness^3 *(TestSetup.width_layered)/ 12
    I_par       = Parylene.thickness^3 *(TestSetup.width_layered )/ 12
    I_steel     = Steel.thickness^3 *(TestSetup.width_steel)/ 12
    
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
        
        k_coh = cohesive_element(CZM.K0, TestSetup.dx,TestSetup.width_layered)
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

function FCZM(CZM0,δ,dN,D)

    A   = CZM0.σ_max/(CZM0.δ_c-CZM0.δ_max)
    dD  = CZM0.m * (δ / CZM0.δ_c) * dN 
    D   = D+dD
    dD= min(D, 1.0)

    K_factor     = (1 - D)
    δ_max_cycle   = (A*CZM0.δ_c)/(A+K_factor*CZM0.K0)
    σ_max_cycle  = δ_max_cycle*K_factor*CZM0.K0
    
    CZM_cycle = CohesiveProperties(σ_max_cycle,δ_max_cycle, CZM0.δ_c, CZM0.m)
    return CZM_cycle,D
end

function fatigue_degradation!(TestSetup,CZM,K, u, cycles, damage, n_dof_steel)
    
   damage_hist      =   zeros(TestSetup.n_elem_layered)
   #δ_last          =   abs(u[n_dof_steel + (TestSetup.n_elem_layered-1)*2 + 1] - u[(TestSetup.start_elem + TestSetup.n_elem_layered- 1-1)*2 + 1])
   a                =   TestSetup.n_elem_layered*TestSetup.dx
   #A               =   CZM.σ_max/(CZM.δ_c-CZM.δ_max)
   #n_coh_element   =   0 # Cohesive element where there is a crack tip
   Gc               =   0
   #δmax_cycle      =   []
   #K_factor_cycle  =   0

   for e in 1:TestSetup.n_elem_layered
        i           = TestSetup.n_elem_layered-e +1 # From last element to the first one
        elem        = TestSetup.start_elem + i - 1
        idx_steel   = (elem-1)*2 + 1
        idx_si      = n_dof_steel + (i-1)*2 + 1
        δ           = abs(u[idx_si] - u[idx_steel])

        if δ > 0
            # Damage per beam element computation
            CZM_cycle,damage[i] = FCZM(CZM,δ,cycles,damage[i])
            #damage[i] += CZM.m * (δ / CZM.δ_c) * cycles
            #damage[i] = min(damage[i], 1.0)

            #push!(damage_hist,CZM.m * (δ_last / CZM.δ_c) * cycles)
            push!(damage_hist,damage[i])
            #K_factor     = (1 - damage[i])
            
            #δmax_cycle   = (A*CZM.δ_c)/(A+K_factor*CZM.K0)
            #σ_max_cycle  = δmax_cycle*K_factor*CZM.K0
            
            σ   =  δ*CZM_cycle.K0

            #if δ>δmax_cycle
               # a=a+TestSetup.dx
               # n_coh_element=n_coh_element+1
            #end

            ## Crack length calculation

            if σ    >   CZM_cycle.σ_max
                a=  a -    TestSetup.dx
                #n_coh_element=i
                Gc   = CZM_cycle.σ_max*CZM_cycle.δ_c/2

            end

             # Update Steifigkeitsmatrix
            idx_steel_full = [(elem-1)*2+1:(elem-1)*2+2; (elem)*2+1:(elem)*2+2]
            idx_si_full = [n_dof_steel + (i-1)*2+1:n_dof_steel + (i-1)*2+2; n_dof_steel + i*2+1:n_dof_steel + i*2+2]

            k_coh_old = cohesive_element(CZM.K0, TestSetup.dx,TestSetup.width_layered)
            k_coh_new = cohesive_element(CZM_cycle.K0, TestSetup.dx,TestSetup.width_layered)
            K[idx_si_full, idx_si_full] += (k_coh_new - k_coh_old)
            K[idx_steel_full, idx_steel_full] += (k_coh_new - k_coh_old)
            K[idx_si_full, idx_steel_full] -= (k_coh_new - k_coh_old)
            K[idx_steel_full, idx_si_full] -= (k_coh_new - k_coh_old)

        end

    end

    #if n_coh_element>0 
        #K_factor_cycle= (1 - damage[n_coh_element])
        #δmax_cycle= 
    #else
        #K_factor_cycle= 1
    #end
   #Gc   = CZM.K0*K_factor_cycle*δmax_cycle*CZM.δ_c/2
    #δmax_cycle= (A*CZM.δ_c)/(A+K_factor_cycle*CZM.K0)
  
    return K,damage_hist,Gc,a
end

function simulate_fatigue(TestSetup,max_force::Float64, n_cycles::Int,Si,Parylene,Steel,CZM)
    damage      = zeros(TestSetup.n_elem_layered)
    Gc_history  = []
    u_history   = Vector{Vector{Float64}}()
    n_dof_steel = length(TestSetup.nodes) * 2  # Defined locally
    
    K_base, F = assemble_system(TestSetup,max_force,Si,Parylene,Steel,CZM)
    #println("Rank of K: ", rank(K_base), " Size: ", size(K_base, 1))
    K = copy(K_base)
    u_initial = K \ F
    #println("Initial Max Steel Deflection (μm): ", maximum(abs.(u_initial[1:2:n_dof_steel])) * 1e6)
    #println("Initial Max Si Deflection (μm): "   , maximum(abs.(u_initial[n_dof_steel+1:2:end])) * 1e6)
    
    damage_history=Vector{Vector{Float64}}()
    a_history=zeros(n_cycles)

    # Computation for each cycle

    for cycle in 1:n_cycles
        u = K \ F   # Displacement(deflection) of beam elements nodes
        push!(u_history, copy(u))
        K = copy(K_base) 

        K,damage_hist, Gc,a = fatigue_degradation!(TestSetup,CZM,K, u, cycle, damage, n_dof_steel)
        push!(damage_history,damage_hist)
        push!(Gc_history,Gc)
        a_history[cycle]=a

        #checkpoints=[1000 2000 3000 4000 5000]
        #if cycle in checkpoints
         #   println("Gc (J/m2): ", Gc)
           # println("a  (mm): ", a*1000)
        #end
        
        
        #if a >  a_history[cycle]
          # a_history[cycle]=a
        #else
           # if cycle>1
             #   a_history[cycle]=a_history[cycle-1]
           # else
             #   a_history[cycle]=0
            #end

        #end

    end
    

    return u_history, damage, damage_history,Gc_history,a_history
end

