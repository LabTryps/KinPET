using DifferentialEquations

function kinetics(model!, parameters, fixed_parameters, c0allExp, tAllExp; exp_batches = Int[], metabolites_to_output = Bool[])
    # Simulate one or more concentration time courses based on a given model and its parameters
    # Supports multiple experiments: just add a new line to c0allExp with the initial concentration of the new experiment and concatenate the desired time points to tAllExp
    #
    # parameters -> parameters[1]=KmA parameters[2]=KmB parameters[3]=KmP parameters[4]=KmQ parameters[5]=Vmax
    # If the experiments were performed with different exp_batches, parameters[6,7,...]=concentration of enzyme/protein, always relative to the first extract    
    # c0allExp -> Initial concentrations (concentrations of the chemical species at time=0). Columns are the chemical species. Each row is one experiment
    # exp_batches -> If you have different exp_batches, this optional vector identifies them by a sequential number, starting at 1
    # t -> Time points at which to calculate concentrations. It should be one column vector concatenating time points for each experiment
    
    C = zeros(length(tAllExp), size(c0allExp,2))
    
    timeCoursesStarts = findall(tAllExp .== 0)  # Assumes that time courses all start at 0
    
    if length(exp_batches) == 0  # There is only one batch
        exp_batches = ones(Int, size(c0allExp, 1))
        param_end = length(parameters)
        parameters = [parameters; 1.0]
    else  # There are other exp_batches, just add the first (reference)
        param_end = length(parameters) - exp_batches[end] + 1
        parameters = [parameters[1:param_end]; 1.0; parameters[param_end+1:end]]
    end

    for i in 1:length(timeCoursesStarts) - 1
        t = tAllExp[timeCoursesStarts[i]:timeCoursesStarts[i + 1] - 1]

        c0 = c0allExp[i, :]
        parametersThisBatch = [parameters[1:param_end]; parameters[param_end + exp_batches[i]]]  # Select the appropriate enzyme concentration
        tspan = (t[1], t[end])
        problem = ODEProblem(model!, c0, tspan, [fixed_parameters; parametersThisBatch])
        solution = solve(problem)
        Cvout = solution(t)
        C[timeCoursesStarts[i]:timeCoursesStarts[i + 1] - 1, :] .= convert(Matrix, Cvout')
    end

    t = tAllExp[timeCoursesStarts[end]:end]

    c0 = c0allExp[end, :]
    parametersThisBatch = [parameters[1:param_end]; parameters[param_end + exp_batches[end]]]  # Select the appropriate enzyme concentration
    tspan = (t[1], t[end])
    problem = ODEProblem(model!, c0, tspan, [fixed_parameters; parametersThisBatch])
    solution = solve(problem)
    Cvout = solution(t)
    C[timeCoursesStarts[end]:end, :] .= convert(Matrix, Cvout')
    
    if length(metabolites_to_output) > 0
        C = C[:, metabolites_to_output]
    end
    C
end
