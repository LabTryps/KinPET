include("kinetics.jl")

using StatsBase, Optim

function fitmodel(model!, fixed_parameters, initial_conc, t, conc, B0, LB, UB; exp_batches = Int[], metabolites_to_output = Bool[], scale = 1.0, options = Optim.Options(), print_solution = false)
    # Find the parameters of a given model that best fit experimental time course data of an enzymatic assay
    #
    # fixed_parameters -> fixed parameters of the model, such as Keq, the equilibrium constant of the reaction
    # t and conc -> time points and corresponding concentrations. Each species in a column. Multiple experiments should be concatenated one below the other. All time courses must start at t = 0
    # If the experiments were performed with different batches, exp_batches should be a sequential number identifying the batch, starting at 1
    # B0 -> Initial guess for the parameters [KmA KmB KmP KmQ Vmax]
    # If the experiments were performed with different batches, B0(6:end) will be the enzyme/protein concentrations relative to the first batch
    # LB and UB: Lower and Upper Bound for the parameters. All must be greater than zero
    # options -> Optimization options for lsqcurvefit
    #
    # Outputs:
    # parameters -> the parameters of the model that best fit the provided data, along with concentration estimates for all exp_batches, relative to that of the first experiment

    # Argument validation
    @assert size(t, 2) == 1
    @assert size(B0, 2) == 1 && all(x -> x > 0, B0)
    @assert size(LB, 2) == 1 && all(x -> x > 0, LB)
    @assert size(UB, 2) == 1 && all(x -> x > 0, UB)
    @assert size(scale) == size(B0) || length(scale) == 1 && all(x -> x > 0, scale)
    @assert size(B0) == size(LB)
    @assert size(B0) == size(UB)

    obj_func(B) = kinetics(model!, B .* scale, fixed_parameters, initial_conc, t, exp_batches=exp_batches, metabolites_to_output=metabolites_to_output)
    loss_func(B) = msd(vec(conc), vec(obj_func(B)))
    solution = optimize(loss_func, LB ./ scale, UB ./ scale, B0 ./ scale, Fminbox(BFGS()), options)
    if print_solution
        println(solution)
    end

    parameters = Optim.minimizer(solution) .* scale
end
