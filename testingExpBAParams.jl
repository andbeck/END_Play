# recreating Binzer 2016 - use case for BEFWM2
#]dev /Users/BOP21HMB/PhD/BEFWM2

# set up
using EcologicalNetworksDynamics # import package
using Plots
using DataFrames
using CSV
using Distributions
using Statistics
#using CairoMakie
using StatsPlots
import Random.seed!


# vectors of variables
T_range = 0:1:40
K_range = 1:1:20

# T and K values for assessment
T_values = 273.15 .+ collect(T_range) # temperature 1-40C
K_int_values = collect(K_range) # intercept of the carrying capacity (eutrophication)

## Food Web Setup ---------------------
n_rep = 1 # number of webs to make (1 for this assessment)

# basal species starting biomass
m0 = 0.01

# constants
websize = 30
ppmr = 100
con = 0.1

# make nrep food web using niche model and Z value 100
seed!(22)

FWs_nb = []

# Food Webs with non-binzer error vals
seed!(22)
for _ in 1:n_rep
    fw_nb = FoodWeb(nichemodel, 30; C=0.1, Z=ppmr, check_disconnected = true, check_cycle = true)
    fw_nb.M *= m0
    push!(FWs_nb, fw_nb)
end

## The testing ---------------

# K vs K_val (the intercept replacement value)

foodweb = FWs_nb[1]
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

    # push data to storage
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

    # push data to storage
     push!(storeT, [T - 273.15, mean(kk), mean(rr), mean(xx), mean(yy)])

 end

storeT

@df storeK StatsPlots.plot(:K, [:meanK, :meanr, :meanx, :meany])
@df storeT StatsPlots.plot(:T, [:meanK, :meanr, :meanx, :meany])

## T diagnotstics
# should be negative expo (check!)
p1 = @df storeT StatsPlots.plot(:T, [:meanK], title = "Neg Expo (Y)", xlabel = "T", ylabel = "K" )

# should be positive expo (check!)
p2 = @df storeT StatsPlots.plot(:T, [:meanr], title = "Pos Expo (Y)", xlabel = "T", ylabel = "r")

# should be negative (FAIL)
p3 = @df storeT StatsPlots.plot(:T, [:meanx], title = "neg 'lin' (FAIL)", xlabel = "T", ylabel = "x")

# should be negative? (FAIL)
p4 = @df storeT StatsPlots.plot(:T, [:meany],title = "neg 'lin' (FAIL)", xlabel = "T", ylabel = "y")

plot(p1,p2,p3,p4, layout =(2,2))

# K diagnotstics

# should be positive linear (check!)
p5 = @df storeK StatsPlots.plot(:K, [:meanK], title = "Pos Linear (Y)", xlabel = "k_intercept", ylabel = "K")

# should be constant (check!)
p6 = @df storeK StatsPlots.plot(:K, [:meanr], title = "Constant (Y)", xlabel = "k_intercept", ylabel = "K")

# should be constant (check!)
p7 = @df storeK StatsPlots.plot(:K, [:meanx], title = "Constant (Y)", xlabel = "k_intercept", ylabel = "K")

# should be constant (check!)
p8 = @df storeK StatsPlots.plot(:K, [:meany], title = "Constant (Y)", xlabel = "k_intercept", ylabel = "K")

plot(p5,p6,p7,p8, layout =(2,2))
