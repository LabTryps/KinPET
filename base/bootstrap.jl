include("fitmodel.jl")

using Random: rand
using LinearAlgebra: dot

function run_bootstrap(bootstrap_replicas, model!, fixed_params, initial_conc, t, conc, B0, LB, UB; exp_batches = Int[], metabolites_to_output = Bool[], scale = 1.0, options = Optim.Options())
    num_of_experiments = size(initial_conc, 1)
    parameters_bootstrap = zeros(Float64, bootstrap_replicas, length(B0))
    exp_choices = zeros(Int, bootstrap_replicas, num_of_experiments)
    
    if length(exp_batches) == 0  # There is only one batch
        exp_batches = ones(Int, size(initial_conc, 1))
        param_end = length(B0)
    else  # There are other exp_batches, just add the first (reference)
        param_end = length(B0) - exp_batches[end] + 1
    end

    Threads.@threads for i in 1:bootstrap_replicas
        # Select experiments for this iteration
        selected_exps = sort(rand(1:num_of_experiments, num_of_experiments))
        exp_choices[i, :] = selected_exps

        # Create the bootstrap sample
        timeCoursesStarts = findall(iszero, t)
        exp_batch_boot = exp_batches[selected_exps]
        t_boot = Vector{Float64}()
        sizehint!(t_boot, length(t))
        conc_boot = Matrix{Float64}(undef,0,size(conc,2))
        initial_conc_boot = initial_conc[selected_exps,:]
        for exp in selected_exps
            if exp < num_of_experiments
                append!(t_boot, t[timeCoursesStarts[exp]:timeCoursesStarts[exp+1]-1])
                conc_boot = [conc_boot; conc[timeCoursesStarts[exp]:timeCoursesStarts[exp+1]-1, :]]
            else
                append!(t_boot, t[timeCoursesStarts[end]:end])
                conc_boot = [conc_boot; conc[timeCoursesStarts[end]:end, :]]
            end
        end

        # Correct exp_batch_boot and remember mapping
        batch_map = Dict(1 => exp_batch_boot[1])
        exp_batch_boot[exp_batch_boot .== exp_batch_boot[1]] .= 1  # If the first sample is not part of batch 1, make it so
        for j in 2:length(exp_batch_boot)
            if exp_batch_boot[j] != exp_batch_boot[j-1]
                batch_map[exp_batch_boot[j-1] + 1] = exp_batch_boot[j]
                exp_batch_boot[exp_batch_boot .== exp_batch_boot[j]] .= exp_batch_boot[j-1] + 1  # Ensure batches follow a sequential order
            end
        end

        # Run the algorithm
        considered_batches = unique(exp_batch_boot)
        if considered_batches[1] == 1
            considered_batches = considered_batches[2:end]
        end
        B0_boot = [B0[1:param_end]; B0[considered_batches .+ (param_end-1)]]
        LB_boot = [LB[1:param_end]; LB[considered_batches .+ (param_end-1)]]
        UB_boot = [UB[1:param_end]; UB[considered_batches .+ (param_end-1)]]
        if length(scale) == 1
            scale_boot = scale
        else
            scale_boot = [scale[1:param_end]; scale[considered_batches .+ (param_end-1)]]
        end
        
        params = fitmodel(model!, fixed_params, initial_conc_boot, t_boot, conc_boot, B0_boot, LB_boot, UB_boot, exp_batches=exp_batch_boot, metabolites_to_output=metabolites_to_output, scale=scale_boot, options=options)

        # Store results in the appropriate row of parameters_bootstrap
        full_params = zeros(Float64, length(B0)) .- 1.0 # Missing batches will be represented with -1.0
        full_params[1:param_end] = params[1:param_end]
        for j in param_end+1:length(params)
            full_params[batch_map[j-param_end]+param_end] = params[j] # Reposition enzime reestimations to original batches
        end
        parameters_bootstrap[i:i, :] .= convert(Matrix, transpose(full_params))
        
        print(string("Finished bootstrap replica ", i, "\n"))
    end
    
    # Vmax and enz amount are linearly dependant, and we take batch 1 as reference
    # For bootstrap replicas where the first batch is not the actual number 1, we have to correct the Vmax and enz amount to match the others
    if exp_batches[end] > 1  # Of course if there is only one batch nothing needs to be done
        # The target will be the mean of the replicas that considered batch 1:
        replicas_with_batch_1 =  parameters_bootstrap[parameters_bootstrap[:,param_end+1] .> 0, param_end+2:end]
        target = zeros(Float64, size(replicas_with_batch_1, 2))
        for i in axes(replicas_with_batch_1)[2]
            target[i] = mean(replicas_with_batch_1[replicas_with_batch_1[:,i] .> 0, i]) # Only calculate using existing values. Missing batches are represented with -1.0
        end

        for i in axes(parameters_bootstrap)[1]
            # If there is no batch 1 in this replica, perform the correction by least squares with the target vector
            if parameters_bootstrap[i,param_end+1] == -1.0
                batch_concentrations = @view parameters_bootstrap[i,param_end+2:end]
                existing_batch_concentrations = batch_concentrations[batch_concentrations .> 0]
                existing_target = target[batch_concentrations .> 0]
                correction_factor = dot(existing_target,existing_batch_concentrations) / dot(existing_batch_concentrations,existing_batch_concentrations)

                parameters_bootstrap[i,param_end] /= correction_factor # Fix Vmax
                parameters_bootstrap[i,param_end+2:end] *= correction_factor # Fix concentrations
                parameters_bootstrap[i,parameters_bootstrap[i,:] .< 0] .= -1.0 # Return non-existing batches to the convention -1
            end
        end
    end

    parameters_bootstrap, exp_choices
end