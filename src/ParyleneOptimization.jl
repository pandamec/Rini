
# Define the predictive function
function StressPrediction(ParyModel,ϵ,dϵdt)
    
    σ_sim = simulate_strain(ParyModel,ϵ,dϵdt)

    return σ_sim
end

# Objective function to minimize (difference between model and experimental data)
function objective_function(ParyModel,σ_exp,ϵ_exp,dϵdt)


    σ_sim    = StressPrediction(ParyModel, ϵ_exp,dϵdt)
    
    # Calculate sum of squared errors
  
    error = sum((σ_sim - σ_exp).^2)

    return error
end

# Main optimization function

function optimize_FCZM(CZM,TestSetup,max_force,N_cycles,fatigueExp; max_iterations=n_it,Si,Parylene,Steel)
    # Initial guess for parameters [Gc0, k1, k2, σ_max]
    
    initial_params = [ParyModel0.k1,ParyModel0.k2,ParyModel0.n]

    lower_bounds = [1e6, 1e6,0]
    upper_bounds = [ 10e6, 10e6, 1e9]

    # Define the objective function with fixed experimental data
    obj(ParyModel_fit) = objective_function(ParyModel_fit,σ_exp,ϵ_exp,dϵdt)
    
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