include("../base/fitmodel.jl")
include("../base/models.jl")
include("../base/bootstrap.jl")

using Plots: plot, plot!

# Validation of the algorithm with simulated data
## Generating data for the fictitious enzyme
# Reaction: A + B <-> P + Q
# Model used: Bi-Bi Reversible Hill Equation (RHE)
# Parameters: [KmA KmB KmP KmQ Vmax] KmQ = Keq*KmA*KmB/KmP = 12
model! = rhe22!
real_params = [20; 40; 200; 6; 2]
Keq = 3

## Time points to be simulated and initial concentrations
# C0: each line is an experiment, with initial conditions [A0 B0 P0 Q0]
# Experimental set starting with no products, saturating B and varying A then vice-versa
# tAllExp: data points for each experiment. All must start with zero
Tsamp = 4; # sampling period
C0 = [2	400	0	0;
      5	400	0	0;
      10	400	0	0;
      20	400	0	0;
      80	400	0	0;
      100	400	0	0;
      150	400	0	0;
      200	400	0	0;
      400   2	0	0;
      400   5	0	0;
      400   10	0	0;
      400   20	0	0;
      400   80	0	0;
      400   100	0	0;
      400   150	0	0;
      400   200	0	0];
tAllExp = repeat(0:Tsamp:1000, outer=16);

# Suppose we repeated the experiment three times, on different days, under the same conditions
C0 = [C0; C0; C0]
tAllExp = [tAllExp; tAllExp; tAllExp]

## Generating the simulated time courses
C = kinetics(model!,real_params,Keq,C0,tAllExp);

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
display(plt)

# In particular, notice that the repeats yield three time courses that overlap perfectly
threeIntervals = repeat(0:Tsamp:1000, outer=3);
experimentsIntervals = [1758:2008; 5774:6024; 9790:10040];

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(threeIntervals, C[experimentsIntervals,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(threeIntervals, C[experimentsIntervals,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(threeIntervals, C[experimentsIntervals,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(threeIntervals, C[experimentsIntervals,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
display(plt)

## Lets try to fit the RHE to these experiments
# Initial guess and lower/upper bounds [KmA KmB KmP Vmax]
# Remember that Keq is a fixed parameter, so will not be optimized
initial_guess = [100;  100;  100;  100;   1  ]
lower_bound   = [1e-6; 1e-6; 1e-6; 1e-6; 1e-6]
upper_bound   = [500;  500;  500;  500;  10  ]

# Options for the algorithm searching for the best fitted parameters
options = Optim.Options(show_trace=false, extended_trace=false, f_tol=1e-5, g_tol=1e-5, time_limit=60, iterations=10^6)
estimated_params = fitmodel(model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)

# Run bootstrap
bootstrap_replicas = 32
parameters_bootstrap, _ = run_bootstrap(bootstrap_replicas, model!, Keq, C0, tAllExp, C, initial_guess, lower_bound, upper_bound, options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)

# Compare the output of the simulations using the optimized parameters with the dataset
sim_C = kinetics(model!,estimated_params,Keq,C0,tAllExp)
plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)

## So far so good. However, enzymes lose activity over time. Let's consider this effect
# Suppose that on day 2, we have 60% of the active enzyme of day 1, and on day 3, 20% of day 1
# First, we'll assign each experiment to a batch according to the day in which it was performed
exp_batches = [fill(1, 16); fill(2, 16); fill(3, 16)]
# Now we'll generate the time courses. The modified enzyme amounts must be declared after the other parameters.
# The amount of the first batch im implicitly 1, and should not be included
real_params = [real_params; 0.6; 0.2]

C = kinetics(model!,real_params,Keq,C0,tAllExp,exp_batches=exp_batches);

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
display(plt)

# Plot the same experiment on different days to make the distinction clearer
plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(threeIntervals, C[experimentsIntervals,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(threeIntervals, C[experimentsIntervals,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(threeIntervals, C[experimentsIntervals,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(threeIntervals, C[experimentsIntervals,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
display(plt)

# If we naively try to use the same strategy as before, ignoring the loss of activity, the parameters will be way off
estimated_params = fitmodel(model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_bootstrap, _ = run_bootstrap(bootstrap_replicas, model!, Keq, C0, tAllExp, C, initial_guess, lower_bound, upper_bound, options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)

# Compare the output of the simulations using the optimized parameters with the dataset
sim_C = kinetics(model!,estimated_params,Keq,C0,tAllExp)

plot!(threeIntervals, sim_C[experimentsIntervals,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Simulated")
plot!(threeIntervals, sim_C[experimentsIntervals,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(threeIntervals, sim_C[experimentsIntervals,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(threeIntervals, sim_C[experimentsIntervals,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
display(plt)

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)

# If we inform the algorythm about the experiment batches, along with estimates for enzyme amounts on them,
# it'll perform the fit accounting for systematic variations in all experiments in the batch
# As before, batch 1 is taken as reference, and its values should not be included in the estimation
initial_guess = [initial_guess; 1;    1  ]
lower_bound   = [lower_bound;  0.01; 0.01]
upper_bound   = [upper_bound;   2;    2  ]

estimated_params = fitmodel(model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,exp_batches=exp_batches,options=options)
parameters_bootstrap, exp_choices = run_bootstrap(bootstrap_replicas, model!, Keq, C0, tAllExp, C, initial_guess, lower_bound, upper_bound, exp_batches=exp_batches, options=options)
# Since the choice of experiments is random for each bootstrap replica, there is always the possibility that no experiments of a given batch were chosen
# If this happens, the corresponding enzyme reestimation for this batch will be marked with -1.0. We can easily exclude these elements by considering only positive values
parameters_boot_mean = zeros(Float64, size(initial_guess))
std_errors = zeros(Float64, size(initial_guess))
for i in axes(parameters_bootstrap)[2]
      current_parameters = parameters_bootstrap[:,i]
      parameters_boot_mean[i] = mean(current_parameters[current_parameters .> 0])
      std_errors[i] = std(current_parameters[current_parameters .> 0])
end
# Estimated params now include a reestimation of the amount of active enzyme, matching the expected values
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)

# Now simulations using the optimized parameters match the original dataset
sim_C = kinetics(model!,estimated_params,Keq,C0,tAllExp,exp_batches=exp_batches)

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(threeIntervals, C[experimentsIntervals,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(threeIntervals, C[experimentsIntervals,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(threeIntervals, C[experimentsIntervals,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(threeIntervals, C[experimentsIntervals,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
plot!(threeIntervals, sim_C[experimentsIntervals,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Simulated")
plot!(threeIntervals, sim_C[experimentsIntervals,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(threeIntervals, sim_C[experimentsIntervals,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(threeIntervals, sim_C[experimentsIntervals,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
display(plt)

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)

# Batches should group together experiments performed on the same day and with the same amount of enzyme
# If two experiments don't match the restriction above, they should be on different batches
