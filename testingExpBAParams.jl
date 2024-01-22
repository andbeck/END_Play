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
using StatsPlots
import Random.seed!


# vectors of variables
T_range = 0:1:40
K_range = 1:1:20

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

## Some testing.

# K vs K_val (the intercept replacement value)

foodweb = FWs_b[1]
Tdefault = 293.15
Kdefault = 3
meanK_K = []
meanK_T = []

storeK = DataFrame(K = [], meanK=[], meanr=[], meanx = [], meany=[])
storeT = DataFrame(T = [], meanK=[], meanr=[], meanx = [], meany=[])


p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2),
                biorates=BioRates(foodweb; d=0))

for (i_K, K) in enumerate(K_int_values)
   set_temperature!(p, Tdefault, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K)))

    # get the values (nothings are placeholders for non plants)
    kk = filter!(x->x!=nothing, p.producer_growth.K)
    rr = filter!(x->x!=0, p.biorates.r)
    xx = filter!(x->x!=0, p.biorates.x)
    yy = filter!(x->x!=0, p.biorates.y)
    

    # average K

    push!(meanK_K,mean(kk))
    push!(storeK, [K, mean(kk), mean(rr), mean(xx), mean(yy)])
end

storeK

# K versus T

for (i_T, T) in enumerate(T_values)
    set_temperature!(p, T, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = Kdefault)))

     # get the values (nothings are placeholders for non plants)
     kk = filter!(x->x!=nothing, p.producer_growth.K)
     rr = filter!(x->x!=0, p.biorates.r)
     xx = filter!(x->x!=0, p.biorates.x)
     yy = filter!(x->x!=0, p.biorates.y)

     # average K

     push!(meanK_T,mean(kk))
     push!(storeT, [T - 273.15, mean(kk), mean(rr), mean(xx), mean(yy)])

 end

storeT

@df storeK StatsPlots.plot(:K, [:meanK, :meanr, :meanx, :meany])
@df storeT StatsPlots.plot(:T, [:meanK, :meanr, :meanx, :meany])

## T diagnotstics
# should be negative expo (check!)
@df storeT StatsPlots.plot(:T, [:meanK])

# should be positive expo (check!)
@df storeT StatsPlots.plot(:T, [:meanr])

# should be negative (FAIL)
@df storeT StatsPlots.plot(:T, [:meanx])

# should be negative? (FAIL)
@df storeT StatsPlots.plot(:T, [:meany])

# K diagnotstics

# should be positive linear (check!)
@df storeK StatsPlots.plot(:K, [:meanK])

# should be constant (check!)
@df storeK StatsPlots.plot(:K, [:meanr])

# should be constant (check!)
@df storeK StatsPlots.plot(:K, [:meanx])

# should be constant (check!)
@df storeK StatsPlots.plot(:K, [:meany])

