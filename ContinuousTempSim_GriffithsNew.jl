using EcologicalNetworksDynamics
using Random, DataFrames, Plots, CSV, Distributions

## OVERVIEW
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
########################################

######################
# Temp options
######################

## reference fixed temperatures
T_fix10 = 10
T_fix40 = 40

## linear increase
# 20 steps from 10c to 40c.
T_lin = collect(range(10, 40 , 20)) # 20 temperature steps

## seasonal flip flops
# need to do this at 10, 25 and 40 with ±1.5 to match T_lin_season
# two temps, each 10x == 20
T_season1 = repeat([10.5, 9.5], 10)
T_season2 = repeat([25.5, 24.5], 10)
T_season3 = repeat([40.5, 39.5], 10)


## sources of variation
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

## test/view distribution of variation
# #extrema - actual range of temperature variations.
# extrema(tvars_superN)
# extrema(tvars_n)
# extrema(tvars_w)

# # visualise distrinbution of temp change; xlims  = extrema of widest
# ptsn = histogram(vec(tvars_superN), xlims = extrema(tvars_w))
# ptn = histogram(vec(tvars_n), xlims = extrema(tvars_w))
# ptw= histogram(vec(tvars_w), xlims = extrema(tvars_w))

# plot(ptsn, ptn, ptw, layout=(3,1))

## Generate Linear + Season
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

## Generate Experimental Design plot
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

##############################################
## MASTER FUNCTION TO DO THE CONTINGENCY WORK
##############################################
simTemp = function(webs, tempSeq)

    # set collection
    df = DataFrame(fw = [], step = [], temp=[], richness = [], stability = [], biomass = [])

    # over replicate (10) webs
    for f in 1:length(webs)

        # start model
        # inital web
        foodweb = webs[f]

        # initial paramters
        p = ModelParameters(foodweb,
                            functional_response=ClassicResponse(foodweb, h=1.2),
                            biorates=BioRates(foodweb; d=0))

        # set initial biomasses
        B0 = zeros(richness(foodweb))
        K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])
        B0[producers(foodweb)] .= K_prod
        B0[1:richness(foodweb) .∉ [producers(foodweb)]] .= K_prod / 8

        # run model across range of T's
        # from specific tempSeq

        for i in 1:size(tempSeq,1)

            #################################################
            ## Organise things for temperature depdendence ##
            #################################################

            # define T for time step i
            TT = tempSeq[i]+273.15 # adjust to Kelvin

            # set K vals
            KK = 1
            K_int = K_int_values[1]
            K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])

            # allocate temperature dependence to rates
            set_temperature!(p, TT, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K_int)))

            ## simulate biomass dynamics for 10 years ##
            out = simulate(p, B0 ,tmax = 3153600000,
                callback = EcologicalNetworksDynamics.ExtinctionCallback(1e-12, p, false),
                adaptive = true,
                dt = 24*60*60,
                saveat = 24*60*60,
            )

            ## collect deets and write to df ##
            rr = richness(out)
            ss = community_cv(out)
            bb = biomass(out).total
            push!(df, [f, i, TT, rr, ss, bb])

            #############################################
            ### identify extinctions in prep for i+1 ####
            #############################################

            # who is extinct
            who_extinct = keys(get_extinct_species(out))

            # a list of 1:n species in current network
            species = 1:size(foodweb.A, 1)
            #extant = setdiff(species, who_extinct)

            # mask extinct species or leave as is if who_extinct is empty
            if !isempty(who_extinct)
                species_idx = setdiff(species, who_extinct) # or findall(x->x ∉ who_extinct, species)
            else
                species_idx = 1:size(foodweb.A, 1)
            end

            ######################################
            ## Update foodweb, B0 and p for t+1 ##
            ######################################

            # subset the matrix
            foodweb = FoodWeb(foodweb.A[species_idx, species_idx])

            # subset and collect the biomass (last value approach vs. mean?)
            # whether the mean or the last value is taken doesn't seem to matter
            B0 = biomass(out, last = 1).species[species_idx]

            # reset the params and bodymass vector with subsetted bodymass
            p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2), biorates=BioRates(foodweb; d=0))
        end
    end
    return(df)
end


################################
## SIMULATIONS USING FUNCTION ##
################################

outLin = simTemp(FWs, T_lin)
outSeason10 = simTemp(FWs, T_season1)
outSeason25 = simTemp(FWs, T_season2)
outSeason40 = simTemp(FWs, T_season3)
outLinSeason = simTemp(FWs, T_lin_season)
out10 = simTemp(FWs, T_fix10)
out40 = simTemp(FWs, T_fix40)

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


# Send to R Visualisation
CSV.write("tempLinRun.csv", outLin)
CSV.write("tempSeason10.csv", outSeason10)
CSV.write("tempSeason25.csv", outSeason25)
CSV.write("tempSeason40.csv", outSeason40)
CSV.write("tempLinSeason.csv", outLinSeason)
CSV.write("temp10Cons.csv", out10)
CSV.write("temp40Cons.csv", out40)

df_LVn = reduce(vcat,df_vector_LV_n)
df_LVlo = reduce(vcat,df_vector_LV_lo)
df_LVhi = reduce(vcat,df_vector_LVS_hi)

CSV.write("tempLinVar_n.csv", df_LVn)
CSV.write("tempLinVar_lo.csv", df_LVlo)
CSV.write("tempLinVar_hi.csv", df_LVhi)

df_LVSn = reduce(vcat,df_vector_LVS_n)
df_LVSlo = reduce(vcat,df_vector_LVS_lo)
df_LVShi = reduce(vcat,df_vector_LVS_hi)

CSV.write("tempLinVarSeason_n.csv", df_LVS_n)
CSV.write("tempLinVarSeason_lo.csv", df_LVS_lo)
CSV.write("tempLinVarSeason_hi.csv", df_LVS_hi)

### Go to exlploreContTempData.R for figures and stats
