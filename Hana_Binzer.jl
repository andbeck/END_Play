# recreating Binzer 2016 - use case for BEFWM2
#]dev /Users/BOP21HMB/PhD/BEFWM2

# set up
using EcologicalNetworksDynamics # import package
using Plots
using DataFrames
using CSV
import Random.seed!

# vectors of variables
T_range = 0:2:40
K_range = [1,10]
n_rep = 10
T_values = 273.15 .+ collect(T_range) # temperature 1-40C
K_int_values = collect(K_range) # intercept of the carrying capacity (eutrophication)
n_T, n_K = length(T_range), length(K_range)

# basal species starting biomass
m0 = 0.01

# make 10 food web using niche model and Z value 100
seed!(22)
FWs = []
for _ in 1:n_rep
    fw = FoodWeb(nichemodel, 30; C=0.1, Z=10)
    fw.M *= m0
    push!(FWs, fw)
end


# dataframe to store results
n_lines = n_T * n_K * n_rep
df = DataFrame(eutrophication=zeros(n_lines), temp=zeros(n_lines), repeat=zeros(n_lines), persistence=zeros(n_lines))

## simulation loop for fig 1b
 # varying intercept of the carrying capacity

Threads.@threads for i_K in 1:n_K
    K_int = K_int_values[i_K]

    # varying temperature
    for (i_T, T) in enumerate(T_values)

        # 100 replicates for each combination of parameters
        for (i_rep, foodweb) in enumerate(FWs)

            println("K_int = $K_int", " T = $T", " rep = $i_rep")


            #generate model parameters
            p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2), biorates=BioRates(foodweb; d=0))

            # set temperature
            set_temperature!(p, T, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K_int)))
           # K = p.environment.K
           # foodweb.M

            # set initial biomasses
            B0 = zeros(EcologicalNetworksDynamics.richness(foodweb))
            K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])
            B0[producers(foodweb)] .= K_prod
            B0[1:EcologicalNetworksDynamics.richness(foodweb) .∉ [producers(foodweb)]] .= K_prod / 8

            # simulate biomass dynamics for 1000 years
            out = simulate(p, B0 , tmax = 31536000000,
                callback = EcologicalNetworksDynamics.ExtinctionCallback(1e-12, p, false),
                adaptive = true,
                dt = 24*60*60,
                saveat = 24*60*60,
                )


            ### Metrics
            # persistence
            persistence = species_persistence(out)

            # each rate being altered by temp


            # push results to dataframe
            i_df = n_T * n_rep * (i_K - 1) + n_rep * (i_T - 1) + i_rep
            df[i_df,:] = [K_int, T, i_rep, persistence]

        end

    end
end

df
CSV.write("Binzer_2016_Z10_coarseK.csv", df)
