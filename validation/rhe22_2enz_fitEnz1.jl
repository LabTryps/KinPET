include("../base/fitmodel.jl")
include("../base/models.jl")
include("../base/bootstrap.jl")

using Plots: plot, plot!
using StatsBase: L2dist
using DSP

# Validation of the algorithm with simulated data
## Generating data for the fictitious enzyme
# Reaction 1: A + B <-> I + Q --> unknown enzyme
# Reaction 2: I + C <-> P + R
# Model used: Bi-Bi Reversible Hill Equation (RHE)
# Fixed parameters: [KmI2 KmP KmR KmC Vmax2 enz2 Keq1 Keq2]
# Parameters to fit to data: [KmA KmB KmI1 KmQ1 Vmax1]
model! = rhe22_2enz_fitEnz1!
real_params = [20 40 200 6 2]';
fixed = [20 200 2 50 2 1 0.5 3]';

## Time points to be simulated and initial concentrations
# C0: each line is an experiment, with initial conditions [A0 B0 P0 Q0]
# Starting with no products, saturating B and varying A
# tAllExp: data points for each experiment. All must start with zero
# We are assuming a single experimental batch, so information about enzyme
# concentration is unnecessary, and will be included in the value of Vmax
Tsamp = 4; # sampling period
C0 = [  2   400   400 0	0 0 0;
        5   400   400 0 0 0 0;
       10	400   400 0 0 0 0;
       20	400   400 0 0 0 0;
       80	400   400 0 0 0 0;
      100	400   400 0 0 0 0;
      150	400   400 0 0 0 0;
      200	400   400 0 0 0 0];
tAllExp = repeat(0:Tsamp:1000, outer=8);

## Generating the simulated time courses
C = kinetics(model!,real_params,fixed,C0,tAllExp);

plt = plot(layout=(2,4), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="C", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="I", subplot=4)
plot!(tAllExp, C[:,5], seriestype=:scatter, ms=2, msw=0, title="P", subplot=5)
plot!(tAllExp, C[:,6], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=6)
plot!(tAllExp, C[:,7], seriestype=:scatter, ms=2, msw=0, title="R", subplot=7);

## Lets try to fit the RHE to these experiments
# Initial guess and lower/upper bounds [KmA KmB KmP Vmax]
# Remember that Keq is a fixed parameter, so will not be optimized
initial_guess = [100;  100;  100;  100;   1  ]
lower_bound   = [1e-6; 1e-6; 1e-6; 1e-6; 1e-6]
upper_bound   = [500;  500;  500;  500;  10  ]

# Assuming that we measure C, there is only information about the
# consumption of C and production of P and R (through stoichiometry)
known_time_courses = [false false true false true false true]';
known_C = C[:,known_time_courses]

# Options for the algorithm searching for the best fitted parameters
options = Optim.Options(show_trace=false, extended_trace=false, f_tol=1e-5, g_tol=1e-5, time_limit=60, iterations=10^6)
estimated_params = fitmodel(model!,fixed,C0,tAllExp,known_C,initial_guess,lower_bound,upper_bound,metabolites_to_output=known_time_courses,options=options)

# Run bootstrap
bootstrap_replicas = 32
parameters_bootstrap, exp_choices = run_bootstrap(bootstrap_replicas, model!, fixed, C0, tAllExp, known_C, initial_guess, lower_bound, upper_bound, metabolites_to_output=known_time_courses, options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)

## Compare the output of the simulations using the optimized parameters with the dataset
sim_C = kinetics(model!,estimated_params,fixed,C0,tAllExp)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
plot!(tAllExp, sim_C[:,5], seriestype=:scatter, ms=2, msw=0, subplot=5)
plot!(tAllExp, sim_C[:,6], seriestype=:scatter, ms=2, msw=0, subplot=6)
plot!(tAllExp, sim_C[:,7], seriestype=:scatter, ms=2, msw=0, subplot=7)
display(plt)

## Complement the dataset with experiments where B and C are varied
C0 = [C0;
      400    2  400 0 0 0 0;
      400    5  400 0 0 0 0;
      400   10	400 0 0 0 0;
      400   20	400 0 0 0 0;
      400   80	400 0 0 0 0;
      400   100	400 0 0 0 0;
      400   150	400 0 0 0 0;
      400   200	400 0 0 0 0;
      400   400  2  0 0 0 0;
      400   400  5  0 0 0 0;
      400   400 10	0 0 0 0;
      400   400 20	0 0 0 0;
      400   400 80	0 0 0 0;
      400   400 100	0 0 0 0;
      400   400 150	0 0 0 0;
      400   400 200	0 0 0 0];
tAllExp = [tAllExp; repeat(0:Tsamp:1000, outer=16)]

C = kinetics(model!,real_params,fixed,C0,tAllExp)

plt = plot(layout=(2,4), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="C", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="I", subplot=4)
plot!(tAllExp, C[:,5], seriestype=:scatter, ms=2, msw=0, title="P", subplot=5)
plot!(tAllExp, C[:,6], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=6)
plot!(tAllExp, C[:,7], seriestype=:scatter, ms=2, msw=0, title="R", subplot=7);

## How well the estimated params match the expanded dataset?
sim_C = kinetics(model!,estimated_params,fixed,C0,tAllExp)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
plot!(tAllExp, sim_C[:,5], seriestype=:scatter, ms=2, msw=0, subplot=5)
plot!(tAllExp, sim_C[:,6], seriestype=:scatter, ms=2, msw=0, subplot=6)
plot!(tAllExp, sim_C[:,7], seriestype=:scatter, ms=2, msw=0, subplot=7)
display(plt)

# Mean squared error between real and simulated data
println(L2dist(C,sim_C))

## Fit the data again, with the expanded dataset and plot the results
known_C = C[:,known_time_courses]
estimated_params = fitmodel(model!,fixed,C0,tAllExp,known_C,initial_guess,lower_bound,upper_bound,metabolites_to_output=known_time_courses,options=options)
parameters_bootstrap, exp_choices = run_bootstrap(bootstrap_replicas, model!, fixed, C0, tAllExp, known_C, initial_guess, lower_bound, upper_bound, metabolites_to_output=known_time_courses, options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)

sim_C = kinetics(model!,estimated_params,fixed,C0,tAllExp)
plt = plot(layout=(2,4), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="C", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="I", subplot=4)
plot!(tAllExp, C[:,5], seriestype=:scatter, ms=2, msw=0, title="P", subplot=5)
plot!(tAllExp, C[:,6], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=6)
plot!(tAllExp, C[:,7], seriestype=:scatter, ms=2, msw=0, title="R", subplot=7)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
plot!(tAllExp, sim_C[:,5], seriestype=:scatter, ms=2, msw=0, subplot=5)
plot!(tAllExp, sim_C[:,6], seriestype=:scatter, ms=2, msw=0, subplot=6)
plot!(tAllExp, sim_C[:,7], seriestype=:scatter, ms=2, msw=0, subplot=7)
display(plt)

# Mean squared error between real and simulated data
println(L2dist(C,sim_C))

###########################################################################
## Repeating the exercise with noisy data
noiseVariance = 10;
C = C + noiseVariance*randn(size(C));

# Assuming we can only measure C, the other curves have to be inferred from stoichiometry
# Reaction 1: A + B <-> I + Q
# Reaction 2: I + C <-> P + R
# C order = [A B C I P Q R]
timeCoursesStarts = findall(tAllExp .== 0)
for experiment_num = 1:length(timeCoursesStarts)-1
    # P(t) = -(C(t) - initialC) + initialP
    C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,5] .= -(C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0[experiment_num,3]) .+ C0[experiment_num,5]
    # R(t) = -(C(t) - initialC) + initialR
    C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,7] .= -(C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0[experiment_num,3]) .+ C0[experiment_num,7]
end
# P(t) = -(C(t) - initialC) + initialP
C[timeCoursesStarts[end]:end,5] .= -(C[timeCoursesStarts[end]:end,3] .- C0[end,3]) .+ C0[end,5];
# R(t) = -(C(t) - initialC) + initialR
C[timeCoursesStarts[end]:end,7] .= -(C[timeCoursesStarts[end]:end,3] .- C0[end,3]) .+ C0[end,7];

plt = plot(layout=(2,4), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="C", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="I", subplot=4)
plot!(tAllExp, C[:,5], seriestype=:scatter, ms=2, msw=0, title="P", subplot=5)
plot!(tAllExp, C[:,6], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=6)
plot!(tAllExp, C[:,7], seriestype=:scatter, ms=2, msw=0, title="R", subplot=7)
display(plt)

## Fitting the model to this extremely noisy dataset
known_C = C[:,known_time_courses]
estimated_params = fitmodel(model!,fixed,C0,tAllExp,known_C,initial_guess,lower_bound,upper_bound,metabolites_to_output=known_time_courses,options=options)
parameters_bootstrap, exp_choices = run_bootstrap(bootstrap_replicas, model!, fixed, C0, tAllExp, known_C, initial_guess, lower_bound, upper_bound, metabolites_to_output=known_time_courses, options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)

sim_C = kinetics(model!,estimated_params,fixed,C0,tAllExp)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
plot!(tAllExp, sim_C[:,5], seriestype=:scatter, ms=2, msw=0, subplot=5)
plot!(tAllExp, sim_C[:,6], seriestype=:scatter, ms=2, msw=0, subplot=6)
plot!(tAllExp, sim_C[:,7], seriestype=:scatter, ms=2, msw=0, subplot=7)
display(plt)

# Mean squared error between sims with estimated parameters and sims with real parameters
sim_C_real = kinetics(model!,real_params,fixed,C0,tAllExp)
println(L2dist(sim_C,sim_C_real))

## Lowpass filter to minimize the noise
# We design a filter of half the desired order because we will use
# filtfilt, which is equivalent to using a filter with double order
Nf = 30;
Fs = 1/Tsamp;
Fpass = 0.001;

lowpass = digitalfilter(Lowpass(Fpass, fs=Fs), FIRWindow(hamming(Nf+1)));

b = DSP.Filters.PolynomialRatio(lowpass, [1.0])
f = range(0, stop=pi, length=1000)
plt = plot(1000*f*Fs/(2*pi), 20*log10.(abs.(freqresp(b,f))), ylims=(-72,Inf), ticks=10, title="Frequency response for lowpass FIR filter", xlabel="Frequency (mHz)", ylabel="Magnitude (dB)", legend=false)
display(plt)

## Apply the filter to the data
C_filt = C
timeCoursesStarts = findall(tAllExp .== 0)
for experiment_num = 1:length(timeCoursesStarts)-1
    # filtfilt runs the filter twice in opposing directions, so there is no delay
    C_filt[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,:] .= filtfilt(lowpass,C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,:])
end
C_filt[timeCoursesStarts[end]:end,:] .= filtfilt(lowpass,C[timeCoursesStarts[end]:end,:]);
plt = plot(layout=(2,4), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C_filt[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C_filt[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C_filt[:,3], seriestype=:scatter, ms=2, msw=0, title="C", subplot=3)
plot!(tAllExp, C_filt[:,4], seriestype=:scatter, ms=2, msw=0, title="I", subplot=4)
plot!(tAllExp, C_filt[:,5], seriestype=:scatter, ms=2, msw=0, title="P", subplot=5)
plot!(tAllExp, C_filt[:,6], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=6)
plot!(tAllExp, C_filt[:,7], seriestype=:scatter, ms=2, msw=0, title="R", subplot=7)

## Run the algorithm on the filtered dataset
known_C = C_filt[:,known_time_courses]
estimated_params = fitmodel(model!,fixed,C0,tAllExp,known_C,initial_guess,lower_bound,upper_bound,metabolites_to_output=known_time_courses,options=options)
parameters_bootstrap, exp_choices = run_bootstrap(bootstrap_replicas, model!, fixed, C0, tAllExp, known_C, initial_guess, lower_bound, upper_bound, metabolites_to_output=known_time_courses, options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)

sim_C = kinetics(model!,estimated_params,fixed,C0,tAllExp)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
plot!(tAllExp, sim_C[:,5], seriestype=:scatter, ms=2, msw=0, subplot=5)
plot!(tAllExp, sim_C[:,6], seriestype=:scatter, ms=2, msw=0, subplot=6)
plot!(tAllExp, sim_C[:,7], seriestype=:scatter, ms=2, msw=0, subplot=7)
display(plt)

## Mean squared errors
# Between sims with estimated parameters and data
println(L2dist(C_filt,sim_C))
# Between sims with real parameters and data
println(L2dist(C_filt,sim_C_real))
# Between sims with estimated parameters and sims with real parameters
println(L2dist(sim_C,sim_C_real))

plt = plot(layout=(2,4), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, sim_C_real[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, sim_C_real[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, sim_C_real[:,3], seriestype=:scatter, ms=2, msw=0, title="C", subplot=3)
plot!(tAllExp, sim_C_real[:,4], seriestype=:scatter, ms=2, msw=0, title="I", subplot=4)
plot!(tAllExp, sim_C_real[:,5], seriestype=:scatter, ms=2, msw=0, title="P", subplot=5)
plot!(tAllExp, sim_C_real[:,6], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=6)
plot!(tAllExp, sim_C_real[:,7], seriestype=:scatter, ms=2, msw=0, title="R", subplot=7)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
plot!(tAllExp, sim_C[:,5], seriestype=:scatter, ms=2, msw=0, subplot=5)
plot!(tAllExp, sim_C[:,6], seriestype=:scatter, ms=2, msw=0, subplot=6)
plot!(tAllExp, sim_C[:,7], seriestype=:scatter, ms=2, msw=0, subplot=7)
display(plt)
