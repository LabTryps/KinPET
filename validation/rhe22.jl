include("../base/fitmodel.jl")
include("../base/models.jl")
include("../base/bootstrap.jl")

using Plots: plot, plot!
using StatsBase: L2dist
using DSP

# Validation of the algorithm with simulated data
## Generating data for the fictitious enzyme
# Reaction: A + B <-> P + Q
# Model used: Bi-Bi Reversible Hill Equation (RHE)
# estimated_params: [KmA KmB KmP KmQ Vf] Vr = Vf*KmP*KmQ/(Keq*KmA*KmB) = 1
model! = rhe22!
real_params = [20; 40; 200; 6; 2]
Keq = 3

## Time points to be simulated and initial concentrations
# C0: each line is an experiment, with initial conditions [A0 B0 P0 Q0]
# Starting with no products, saturating B and varying A
# tAllExp: data points for each experiment. All must start with zero
# We are assuming a single experimental batch, so information about enzyme
# concentration is unnecessary, and will be included in the value of Vmax
Tsamp = 4; # sampling period
C0 = [2	    400	0	0;
      5	    400	0	0;
      10	400	0	0;
      20	400	0	0;
      80	400	0	0;
      100	400	0	0;
      150	400	0	0;
      200	400	0	0];
tAllExp = repeat(0:Tsamp:1000, outer=8);

## Generating the simulated time courses
C = kinetics(model!,real_params,Keq,C0,tAllExp);

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4);

## Lets try to fit the RHE to these experiments
# Initial guess and lower/upper bounds [KmA KmB KmP Vmax]
# Remember that Keq is a fixed parameter, so will not be optimized
initial_guess = [100;  100;  100;  100;   1  ]
lower_bound   = [1e-6; 1e-6; 1e-6; 1e-6; 1e-6]
upper_bound   = [500;  500;  500;  500;  10  ]

# Options for the algorithm searching for the best fitted estimated_params
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
# Calculating Vr
println(estimated_params[3]*estimated_params[4]*estimated_params[5]/(Keq*estimated_params[1]*estimated_params[2]))

## Compare the output of the simulations using the optimized estimated_params with the dataset
sim_C = kinetics(model!,estimated_params,Keq,C0,tAllExp)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)

## Complement the dataset with experiments where B is varied
C0 = [C0;
        400 2	0	0;
        400 5	0	0;
        400 10	0	0;
        400 20	0	0;
        400 80	0	0;
        400 100	0	0;
        400 150	0	0;
        400 200	0	0];
tAllExp = [tAllExp; repeat(0:Tsamp:1000, outer=8)]

C = kinetics(model!,real_params,Keq,C0,tAllExp)

## Fit the data again, with the expanded dataset and plot the results
estimated_params = fitmodel(model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_bootstrap, _ = run_bootstrap(bootstrap_replicas,model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)
# Calculating Vr
println(estimated_params[3]*estimated_params[4]*estimated_params[5]/(Keq*estimated_params[1]*estimated_params[2]))

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

###########################################################################
## Repeating the exercise with noisy data
# First, consider the experiments varying the initial concentration of A
# Try out different values for the variance and see how they affect the analysis
noiseVariance = 10;
C0 = [2	    400	0	0;
      5	    400	0	0;
      10	400	0	0;
      20	400	0	0;
      80	400	0	0;
      100	400	0	0;
      150	400	0	0;
      200	400	0	0];
tAllExp = repeat(0:Tsamp:1000, outer=8);
C = kinetics(model!,real_params,Keq,C0,tAllExp);
C = C + noiseVariance*randn(size(C));

# Assuming we can only measure P, the other curves have to be inferred from stoichiometry
timeCoursesStarts = findall(tAllExp .== 0)
for experiment_num = 1:length(timeCoursesStarts)-1
    # Q(t) = P(t) - initialP + initialQ
    C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,4] .= C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0[experiment_num,3] .+ C0[experiment_num,4]
    # A(t) = -(P(t) - initialP) + initialA
    C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,1] .= -(C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0[experiment_num,3]) .+ C0[experiment_num,1]
    # B(t) = -(P(t) - initialP) + initialB
    C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,2] .= -(C[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0[experiment_num,3]) .+ C0[experiment_num,2]
end
# Q(t) = P(t) - initialP + initialQ
C[timeCoursesStarts[end]:end,4] .= C[timeCoursesStarts[end]:end,3] .- C0[end,3] .+ C0[end,4];
# A(t) = -(P(t) - initialP) + initialA
C[timeCoursesStarts[end]:end,1] .= -(C[timeCoursesStarts[end]:end,3] .- C0[end,3]) .+ C0[end,1];
# B(t) = -(P(t) - initialP) + initialB
C[timeCoursesStarts[end]:end,2] .= -(C[timeCoursesStarts[end]:end,3] .- C0[end,3]) .+ C0[end,2];

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
display(plt)

## Optimize the estimated_params to fit the data
# This time we can make use of the confidence interval estimation
estimated_params = fitmodel(model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_bootstrap, _ = run_bootstrap(bootstrap_replicas,model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)
# Calculating Vr
println(estimated_params[3]*estimated_params[4]*estimated_params[5]/(Keq*estimated_params[1]*estimated_params[2]))

sim_C = kinetics(model!,estimated_params,Keq,C0,tAllExp)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)

## Including the extra experiments, where A is satutrating and B is varied
C0_complement = [400 2	    0	0;
                 400 5	    0	0;
                 400 10	    0	0;
                 400 20	    0	0;
                 400 80	    0	0;
                 400 100	0	0;
                 400 150	0	0;
                 400 200	0	0];
t_complement = repeat(0:Tsamp:1000, outer=8)
C_complement = kinetics(model!,real_params,Keq,C0_complement,t_complement);
C_complement = C_complement + noiseVariance*randn(size(C_complement));

# Assuming we can only measure P, the other curves have to be inferred from stoichiometry
timeCoursesStarts = findall(tAllExp .== 0)
for experiment_num = 1:length(timeCoursesStarts)-1
    # Q(t) = P(t) - initialP + initialQ
    C_complement[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,4] .= C_complement[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0_complement[experiment_num,3] .+ C0_complement[experiment_num,4]
    # A(t) = -(P(t) - initialP) + initialA
    C_complement[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,1] .= -(C_complement[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0_complement[experiment_num,3]) .+ C0_complement[experiment_num,1]
    # B(t) = -(P(t) - initialP) + initialB
    C_complement[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,2] .= -(C_complement[timeCoursesStarts[experiment_num]:timeCoursesStarts[experiment_num+1]-1,3] .- C0_complement[experiment_num,3]) .+ C0_complement[experiment_num,2]
end
# Q(t) = P(t) - initialP + initialQ
C_complement[timeCoursesStarts[end]:end,4] .= C_complement[timeCoursesStarts[end]:end,3] .- C0_complement[end,3] .+ C0_complement[end,4];
# A(t) = -(P(t) - initialP) + initialA
C_complement[timeCoursesStarts[end]:end,1] .= -(C_complement[timeCoursesStarts[end]:end,3] .- C0_complement[end,3]) .+ C0_complement[end,1];
# B(t) = -(P(t) - initialP) + initialB
C_complement[timeCoursesStarts[end]:end,2] .= -(C_complement[timeCoursesStarts[end]:end,3] .- C0_complement[end,3]) .+ C0_complement[end,2];

C0 = [C0; C0_complement];
C = [C; C_complement];
tAllExp = [tAllExp; t_complement];

## Fit the data again, with the expanded dataset and plot the resultsplt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
estimated_params = fitmodel(model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_bootstrap, _ = run_bootstrap(bootstrap_replicas,model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)
# Calculating Vr
println(estimated_params[3]*estimated_params[4]*estimated_params[5]/(Keq*estimated_params[1]*estimated_params[2]))

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

## Calculate the mean squared error
sim_C_real = kinetics(model!,real_params,Keq,C0,tAllExp);

# Between sims with estimated estimated_params and data
println(L2dist(C,sim_C))
# Between sims with real estimated_params and data
println(L2dist(C,sim_C_real))
# Between sims with estimated estimated_params and sims with real estimated_params
println(L2dist(sim_C,sim_C_real))

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, sim_C_real[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, sim_C_real[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, sim_C_real[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, sim_C_real[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)

## As an attempt to mitigate the noise, we'll design a lowpass filter
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
plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, C_filt[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, C_filt[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, C_filt[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, C_filt[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)

## Run the algorithm on the filtered dataset
estimated_params = fitmodel(model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_bootstrap, _ = run_bootstrap(bootstrap_replicas,model!,Keq,C0,tAllExp,C,initial_guess,lower_bound,upper_bound,options=options)
parameters_boot_mean = mean(parameters_bootstrap, dims=1)
std_errors = std(parameters_bootstrap, dims=1)
println(estimated_params)
println(parameters_boot_mean)
println(std_errors)
# Calculating Vr
println(estimated_params[3]*estimated_params[4]*estimated_params[5]/(Keq*estimated_params[1]*estimated_params[2]))

sim_C = kinetics(model!,estimated_params,Keq,C0,tAllExp)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)

## Calculate the mean squared error
# Between sims with estimated estimated_params and data
println(L2dist(C_filt,sim_C))
# Between sims with real estimated_params and data
println(L2dist(C_filt,sim_C_real))
# Between sims with estimated estimated_params and sims with real estimated_params
println(L2dist(sim_C,sim_C_real))

plt = plot(layout=(2,2), xlabel="t", ylabel="Concentration", legend=false)
plot!(tAllExp, sim_C_real[:,1], seriestype=:scatter, ms=2, msw=0, title="A", subplot=1, legend=true, labels="Real")
plot!(tAllExp, sim_C_real[:,2], seriestype=:scatter, ms=2, msw=0, title="B", subplot=2)
plot!(tAllExp, sim_C_real[:,3], seriestype=:scatter, ms=2, msw=0, title="P", subplot=3)
plot!(tAllExp, sim_C_real[:,4], seriestype=:scatter, ms=2, msw=0, title="Q", subplot=4)
plot!(tAllExp, sim_C[:,1], seriestype=:scatter, ms=2, msw=0, subplot=1, legend=true, labels="Simulated")
plot!(tAllExp, sim_C[:,2], seriestype=:scatter, ms=2, msw=0, subplot=2)
plot!(tAllExp, sim_C[:,3], seriestype=:scatter, ms=2, msw=0, subplot=3)
plot!(tAllExp, sim_C[:,4], seriestype=:scatter, ms=2, msw=0, subplot=4)
display(plt)
