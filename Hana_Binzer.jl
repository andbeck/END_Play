# recreating Binzer 2016 - use case for BEFWM2
#]dev /Users/BOP21HMB/PhD/BEFWM2

# set up
using EcologicalNetworksDynamics # import package
using Plots
using DataFrames
using CSV
using Distributions
using Statistics
using CairoMakie
import Random.seed!


# vectors of variables
T_range = 0:5:40
K_range = 1:10:20

n_rep = 2

T_values = 273.15 .+ collect(T_range) # temperature 1-40C
K_int_values = collect(K_range) # intercept of the carrying capacity (eutrophication)
n_T, n_K = length(T_range), length(K_range)

# basal species starting biomass
m0 = 0.01

# constants
websize = 30
ppmr = 100
con = 0.1

# make nrep food web using niche model and Z value 100
seed!(22)

# define variation in species Masses
mu = 0
sd = 0.01
d = Normal(mu, sd)
randM = rand(d, websize)

FWs_nb = []
FWs_b = []

# Food Webs with non-binzer error vals
seed!(22)
for _ in 1:n_rep
    #fw_nb = FoodWeb(nichemodel, 30; C=0.1, Z=ppmr)
    fw_nb = FoodWeb(nichemodel, 30; C=0.1, Z=ppmr, check_disconnected = true, check_cycle = true)
    fw_nb.M *= m0
    push!(FWs_nb, fw_nb)
end

# Food Webs with Binzer error vals
seed!(22)
for _ in 1:n_rep
    #fw_b = FoodWeb(nichemodel, 30; C=0.1, Z=ppmr)
    fw_b = FoodWeb(nichemodel, 30; C=0.1, Z=ppmr, check_disconnected = true, check_cycle = true)
    fw_b.M = m0 * 100 .^ (trophic_levels(fw_b) .- 1 + randM)
    push!(FWs_b, fw_b)
end


# test Mass variation
# species in same trophic levels don't have == mass
DataFrame(M_nb = FWs_nb[1].M,  TL_nb = trophic_levels(FWs_nb[1]),
            M_binz = FWs_b[1].M, TL_binz = trophic_levels(FWs_b[1]))

# dataframe to store results
n_lines = n_T * n_K * n_rep
df = DataFrame(eutrophication=zeros(n_lines), temp=zeros(n_lines), repeat=zeros(n_lines), persistence=zeros(n_lines))

# #= ## Some testing.

# # K vs K_val (the intercept replacement value)

# foodweb = FWs_b[1]
# Tdefault = 293.15
# Kdefault = 3
# meanK_K = []
# meanK_T = []

# p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2),
#                 biorates=BioRates(foodweb; d=0))

# for (i_K, K) in enumerate(K_int_values)
#    set_temperature!(p, Tdefault, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K)))

#     # get the values (nothings are placeholders for non plants)
#     kk=filter!(x->x!=nothing, p.producer_growth.K)

#     # average K

#     push!(meanK_K,mean(kk))
# end

# # K versus T

# for (i_T, T) in enumerate(T_values)
#     set_temperature!(p, T, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = Kdefault)))

#      # get the values (nothings are placeholders for non plants)
#      kk=filter!(x->x!=nothing, p.producer_growth.K)
#      # average K

#      push!(meanK_T,mean(kk))
#  end


# p1 = plot(1:20, meanK_K)
# p2 = plot(T_values .- 273.15, meanK_T)
# plot(p1, p2)

#
## simulation loop for fig 1b
 # varying intercept of the carrying capacity

Threads.@threads for i_K in 1:n_K
    K_int = K_int_values[i_K]
    # varying temperature
    for (i_T, T) in enumerate(T_values)

        # 100 replicates for each combination of parameters
        for (i_rep, foodweb) in enumerate(FWs_b)

            println("K_int = $K_int", " T = $T", " rep = $i_rep")

            #generate model parameters
            # sets K = 1 for all producers
            p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2),
                biorates=BioRates(foodweb; d=0))

            # set temperature and intercept of K according to Binzer
            # "eutrophication (varying the intercept of the carrying capacity, d, between 1 and 20 in steps of 1
            # to manipulate the energy input into the food webs)
            # applies intecept and exponents for K, r, a, h and x
            set_temperature!(p, T, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K_int)))

            # set initial biomasses
            B0 = zeros(EcologicalNetworksDynamics.richness(foodweb))
            K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])

            # debugging broadcasting error.
            println(K_prod)

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

#Makie plotting

gdf = groupby(df, [:eutrophication, :temp])
meanPers = combine(gdf, :persistence => mean)

fig = Figure()
ax = Axis(fig[1, 1])
pp = CairoMakie.heatmap!(ax, meanPers[!,:eutrophication], meanPers[!,:temp], meanPers[!,:persistence_mean])
Colorbar(fig[1,2], pp)
fig


# write data for R plotting
CSV.write("Binzer_2016_Z100_ExpBASimple.csv", df)
