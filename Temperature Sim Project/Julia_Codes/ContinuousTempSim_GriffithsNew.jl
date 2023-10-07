###############################################################
## OVERVIEW
###############################################################

# This script implements work to evaluate the impacts of continuous increases in temperature on the biomass, stablility
# and diversity of ecological communities.  The UNIQUE aspect of this work is evaluating how a single change in temperature generates
# changes in biodviersity, biomass, stability and structure of a community that represents a contingecy
# for the next change: the impacts of changes at t+2 are contingent on changes that occured between t and t+1.

## We simulate multiple scenarios

## REFERENCE POINTS
# Fixed Temperatures at 10 and 40 (endpoints)

## CONTINGENCY SCENARIOS
# Seasonal like variation at 10, 25 and 40 (constant mean, cycling temps)
# Linear Increase from 10 to 40 (20 steps of temperature change)
# # ToDo: is this a sensible range?
# # original was 3 x 4 degree changes to mimic IPCC.
# Linear Increase with Seasonal like cycles
# Linear Increase with stochastic variation
# # ±5C max, random, via truncated t distribution
# # ±10C max,  random, via truncated t distribution
# # ± 15C max, random, via truncated t distribution

## Method
# Implement END model starting at 10C, for temp x and run to eq
# collect biomass, bidiversity and stability data
# update network with extinction events
# Implment END model with updated network at next temp
# Repeat

## NOTES ###############################
# What to do about disconnected species
# Better to do 4c changes from 12, 18, 25c?
########################################

######################
# Packages
######################

using EcologicalNetworksDynamics
using Random, DataFrames, Plots, CSV, Distributions

###############################################################
## Include Function to loop over webs and temperature values ##
## Imports SimTemp function                                  ##
###############################################################

include("TemperatureWebLoopFunction.jl")


######################
# Temp options
######################

## reference fixed temperatures
T_fix10 = 10
T_fix40 = 40

## linear increase
# 20 steps from 10c to 40c.
T_lin = collect(range(10, 40 , 20)) # 20 temperature steps

T_lin_10_14 = collect(range(10, 14, 20))

## seasonal flip flops
# need to do this at 10, 25 and 40 with ±1.5 to match T_lin_season
# two temps, each 10x == 20
T_season1 = repeat([10.5, 9.5], 10)
T_season2 = repeat([25.5, 24.5], 10)
T_season3 = repeat([40.5, 39.5], 10)

#########################
## sources of variation
#########################

# using t versus normal to change frequency of things in tails
Random.seed!(8675309)
# truncated normal with max change of ±5
# nn = truncated(Normal(0,2), -5,5)
# trucated t super narrow ±5 change
tt_supernarrow = truncated(TDist(2), -5, 5)
# low levels extreme
# truncated t with max change of ±10
tt_narrow = truncated(TDist(2), -10, 10)
# high levels extreme
# truncated t with max change of ±10
tt_wide = truncated(TDist(2), -20, 20)

# generate 20 random values x 50 replicates to add to linear and season
tvars_superN = rand(tt_supernarrow, length(T_lin), 50)
tvars_n = rand(tt_narrow, length(T_lin), 50)
tvars_w= rand(tt_wide, length(T_lin), 50)

# test/view distribution of variation
#extrema - actual range of temperature variations.
extrema(tvars_superN)
extrema(tvars_n)
extrema(tvars_w)

# visualise distrinbution of temp change; xlims  = extrema of widest
ptsn = histogram(vec(tvars_superN), xlims = extrema(tvars_w))
ptn = histogram(vec(tvars_n), xlims = extrema(tvars_w))
ptw= histogram(vec(tvars_w), xlims = extrema(tvars_w))

# distributions of temp variation
plot(ptsn, ptn, ptw, layout=(3,1))

###################################
## Generate Linear + Cycle (Season)
###################################

# ensure length is correct
ll = trunc(Int, length(T_lin)/2)
# generate cycle on linear increase
# add 1.5 to start, than -3: generates 0.5c changes on line
var_season = repeat([1.5, -3], ll)
# add seasonal to T
T_lin_season = collect(T_lin) + var_season

## Generate Linear + Variation
# add randoms to T
T_lin_varNorm = T_lin .+ tvars_superN # super narrow
T_lin_varloExt = T_lin .+ tvars_n # narrow
T_lin_varhighExt = T_lin .+ tvars_w # wide

## Linear + Season + Variation
T_lin_season_varNorm = T_lin_season .+ tvars_superN
T_lin_season_varloExt = T_lin_season .+ tvars_n
T_lin_season_varhighExt= T_lin_season .+ tvars_w

#######################################
## Generate Experimental Design plot
#######################################

p1 = plot(1:1:20, T_lin, legend = false,
    title = "L")
p2 = plot(1:1:20, [T_season1, T_season2, T_season3],
    title = "S", legend = false)
p3 = plot(1:1:20, T_lin_season, legend = false,
    title = "L+S")
# linear with variations
p4 = plot(1:1:20, T_lin_varNorm, legend = false,
    title = "L+V superLoExt", ylims = extrema(T_lin_season_varhighExt))

p5 = plot(1:1:20, T_lin_varloExt, legend = false,
    title = "L+V loExt", ylims = extrema(T_lin_season_varhighExt))

p6 = plot(1:1:20, T_lin_varhighExt, legend = false,
    title = "L+V hiExt", ylims = extrema(T_lin_season_varhighExt))

# linear with season and variations
p7 = plot(1:1:20, T_lin_season_varNorm, legend = false,
    title = "L+V+S superLoExt", ylims = extrema(T_lin_season_varhighExt))

p8 = plot(1:1:20, T_lin_season_varloExt, legend = false,
    title = "L+V+S loExt", ylims = extrema(T_lin_season_varhighExt))

p9 = plot(1:1:20, T_lin_season_varhighExt, legend = false,
    title = "L+V+S hiExt", ylims = extrema(T_lin_season_varhighExt))

plot(p1, p2, p3, 
    p4, p5, p6,
    p7, p8, p9, layout = (3,3))


######################
# K Stuff
######################

# Currently only using K = 1

K_range = 1.:1:20 # K valus

# organise some K thing.s
K_int_values = collect(K_range) # intercept of the carrying capacity (eutrophication)

# basal species starting biomass
m0 = 0.01

######################
# food webs from niche model
######################

# make 25 food web using niche model and Z value 100
# species richness of 30
# connetance  = 0.1

Random.seed!(22)
web_rep = 25 # reps
FWs = []
for _ in 1:web_rep
    fw = FoodWeb(nichemodel, 30; C=0.1, Z=100)
    fw.M *= m0
    push!(FWs, fw)
end

# n vals for counting
n_T, n_K = 20, 20
n_web = length(FWs)

# dataframe to store results
n_lines = n_web * n_T # 20 temps, 25 webs

################################
## SIMULATIONS USING FUNCTION ##
################################

# reference
out10 = simTemp(FWs, T_fix10)
out40 = simTemp(FWs, T_fix40)

# season at three temps
outSeason10 = simTemp(FWs, T_season1)
outSeason25 = simTemp(FWs, T_season2)
outSeason40 = simTemp(FWs, T_season3)

# linear
outLin = simTemp(FWs, T_lin)
outLin_10_14 = simTemp(FWs, T_lin_10_14)

# linear with seasons
outLinSeason = simTemp(FWs, T_lin_season)

## Linear with Variation (reps)

# super narrow
df_vector_LV_n = Any[]
for i in 1:size(T_lin_varNorm, 2)
    push!(df_vector_LV_n, simTemp(FWs, T_lin_varNorm[:,i]))
end

#wider variation
df_vector_LV_lo = Any[]
for i in 1:size(T_lin_varloExt, 2)
    push!(df_vector_LV_lo, simTemp(FWs, T_lin_varloExt[:,i]))
end

# extremes
df_vector_LV_hi = Any[]
for i in 1:size(T_lin_varhighExt, 2)
    push!(df_vector_LV_hi, simTemp(FWs, T_lin_varhighExt[:,i]))
end

##Linear with Season and Variation (reps)

# narrow
df_vector_LVS_n = Any[]
for i in 1:size(T_lin_season_varNorm, 2)
    push!(df_vector_LVS_n, simTemp(FWs, T_lin_season_varNorm[:,i]))
end

# wider variation
df_vector_LVS_lo = Any[]
for i in 1:size(T_lin_season_varloExt, 2)
    push!(df_vector_LVS_lo, simTemp(FWs, T_lin_season_varloExt[:,i]))
end

#extremes
df_vector_LVS_hi = Any[]
for i in 1:size(T_lin_season_varhighExt, 2)
    push!(df_vector_LVS_hi, simTemp(FWs, T_lin_season_varhighExt[:,i]))
end


## Send to R Visualisation

# reference
CSV.write("Data4R/temp10Cons.csv", out10)
CSV.write("Data4R/temp40Cons.csv", out40)

# seasons
CSV.write("Data4R/tempSeason10.csv", outSeason10)
CSV.write("Data4R/tempSeason25.csv", outSeason25)
CSV.write("Data4R/tempSeason40.csv", outSeason40)

# linear
CSV.write("Data4R/tempLinRun.csv", outLin)

# linear with season
CSV.write("Data4R/tempLinSeason.csv", outLinSeason)

# linear with variation
# requires stacking matrix of 50 random changes using reduce and vcat
df_LVn = reduce(vcat,df_vector_LV_n)
df_LVlo = reduce(vcat,df_vector_LV_lo)
df_LVhi = reduce(vcat,df_vector_LV_hi)

CSV.write("Data4R/tempLinVar_n.csv", df_LVn)
CSV.write("Data4R/tempLinVar_lo.csv", df_LVlo)
CSV.write("Data4R/tempLinVar_hi.csv", df_LVhi)

# linear with season and variation
# requires stacking matrix of 50 random changes using reduce and vcat
df_LVSn = reduce(vcat,df_vector_LVS_n)
df_LVSlo = reduce(vcat,df_vector_LVS_lo)
df_LVShi = reduce(vcat,df_vector_LVS_hi)

CSV.write("Data4R/tempLinVarSeason_n.csv", df_LVSn)
CSV.write("Data4R/tempLinVarSeason_lo.csv", df_LVSlo)
CSV.write("Data4R/tempLinVarSeason_hi.csv", df_LVShi)

###########################################################
### Go to exlploreContTempData.R for figures and stats ####
###########################################################

# @df outLin_10_14 plot(:step, :richness, group = :fw)