
# Define the cohesive zone model function
function FCZMPrediction(CZM_fit, TestSetup,n_cycles,max_force,Si,Parylene,Steel)
    CZM = CohesiveProperties(CZM_fit[1],CZM_fit[2],CZM_fit[3],CZM_fit[4])
    u_hist, damage,damage_history,Gc_history,a_history = simulate_fatigue(TestSetup,max_force, n_cycles,Si,Parylene,Steel,CZM)

    crack_length=a_history[end]

    n_dof_steel = length(TestSetup.nodes) * 2
    u_last = u_hist[end]
    steel_deflection = u_last[1:2:n_dof_steel] 

        si_deflection_last = u_last[(TestSetup.start_elem+TestSetup.n_elem_layered)*2 + 1]  - steel_deflection[TestSetup.start_elem]
        parylene_deflection_last = u_last[n_dof_steel + (TestSetup.n_elem_layered)*2 + 1]
        s=abs(si_deflection_last - parylene_deflection_last)

    fatigueModel=FatigueData(crack_length,s)
    println("CZM(K0,σ_max,δ_max,δ_c,m,G_c): ", CZM)
    return fatigueModel
end

# Objective function to minimize (difference between model and experimental data)
function objective_function(CZM_fit, TestSetup,max_force,n_cycles,fatigueExp,Si,Parylene,Steel)

    fatigueModel    = FCZMPrediction(CZM_fit, TestSetup,n_cycles,max_force,Si,Parylene,Steel)
    
    # Calculate sum of squared errors
  
    error = (fatigueModel.crack_length- fatigueExp.crack_length)^2 + (fatigueModel.separation- fatigueExp.separation)^2 
    println("error: ", error)
    return error
end

# Main optimization function

function optimize_FCZM(CZM,TestSetup,max_force,N_cycles,fatigueExp; max_iterations=n_it,Si,Parylene,Steel)
    # Initial guess for parameters [Gc0, k1, k2, σ_max]
    
    initial_params = [CZM.σ_max, CZM.δ_max,CZM.δ_c, CZM.m]
    
    # Bounds for parameters (adjust as needed)
    lower_bounds = [0.1, 0.0001,  0.0001,1e-10]
    upper_bounds = [ 2e6, 0.001, 0.001, 1e-5]
    
    # Define the objective function with fixed experimental data
    obj(CZM_fit) = objective_function(CZM_fit, TestSetup,max_force,N_cycles,fatigueExp,Si,Parylene,Steel)
    
    # Perform optimization using Nelder-Mead or L-BFGS (you can switch methods)
    result = optimize(obj, 
                     lower_bounds, 
                     upper_bounds, 
                     initial_params, 
                     Fminbox(LBFGS()),
                     Optim.Options(iterations=max_iterations, show_trace=true))
    
    # Extract optimized parameters
    optimized_params = Optim.minimizer(result)
    minimum_error = Optim.minimum(result)
    
    return optimized_params, minimum_error
end