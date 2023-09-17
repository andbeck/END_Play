using EcologicalNetworksDynamics
using Random, DataFrames, Plots, CSV, Distributions

######################
# Temp options
######################

# min and max
T_fix20 = 20
T_fix40 = 40

# linear increase
T_lin = collect(range(20, 40 , 20)) # 20 temperature steps

# seasonal flip flops
# need to do this at 20, 30 and 40 with ±1.5 to match T_lin_season
T_season1 = repeat([19.5, 20.5], 10)
T_season2 = repeat([30.5, 29.5], 10)
T_season3 = repeat([40.5, 39.5], 10)


# linear with variation
# consider t versus normal to change frequency of things in tails
Random.seed!(123)
d = Normal(1.5,3)

# 10 rvars
rvars = rand(d, length(T_lin), 10)
T_lin_var = T_lin .+ rvars

# linear with 'season'
ll = trunc(Int, length(T_lin)/2)
# ensure length is correct
var_season = repeat([1.5, -3], ll)

T_lin_season = collect(T_lin) + var_season

T_lin_season_var = T_lin_season .+ rvars

p1 = plot(1:1:20, T_lin, legend = false,
    title = "L")
p2 = plot(1:1:20, [T_season1, T_season2, T_season3],
    title = "S", legend = false)
p3 = plot(1:1:20, T_lin_season, legend = false,
    title = "L+S")
p4 = plot(1:1:20, T_lin_var, legend = false,
    title = "L+V")
p5 = plot(1:1:20, T_lin_season_var, legend = false,
    title = "L+S+V")

plot(p1, p2, p3, p4, p5, layout = (2,3))

######################
# K Stuff
######################

K_range = 1.:1:20 # K valus

# organise some K thing.s
K_int_values = collect(K_range) # intercept of the carrying capacity (eutrophication)

# basal species starting biomass
m0 = 0.01

######################
# food webs from niche model
######################

# make 10 food web using niche model and Z value 100
Random.seed!(22)
web_rep = 10 # reps
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
n_lines = n_web * n_T # 20 temps, 10 webs

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
outSeason20 = simTemp(FWs, T_season1)
outSeason30 = simTemp(FWs, T_season2)
outSeason40 = simTemp(FWs, T_season3)
outLinSeason = simTemp(FWs, T_lin_season)
out20 = simTemp(FWs, T_fix20)
out40 = simTemp(FWs, T_fix40)

#Linear with Variation (reps)
df_vector_LV = Any[]
for i in 1:size(T_lin_var, 2)
    push!(df_vector_LV, simTemp(FWs, T_lin_var[:,i]))
end


#Linear with Season and Variation (reps)
df_vector_LVS = Any[]
for i in 1:size(T_lin_season_var, 2)
    push!(df_vector_LVS, simTemp(FWs, T_lin_season_var[:,i]))
end

# Send to R Visualisation
CSV.write("tempLinRun.csv", outLin)
CSV.write("tempSeason20.csv", outSeason20)
CSV.write("tempSeason30.csv", outSeason30)
CSV.write("tempSeason40.csv", outSeason40)
CSV.write("tempLinSeason.csv", outLinSeason)
CSV.write("temp20Cons", out20)
CSV.write("temp40Cons", out40)
df_LV = reduce(vcat,df_vector_LV)
CSV.write("tempLinVar.csv", df_LV)
df_LV = reduce(vcat,df_vector_LVS)
CSV.write("tempLinVarSeason.csv", df_LVS)


# p1 = plot(df[!,:temp] .-273.15, df[!,:richness])
# ylabel!("Richness")
# p2 = plot(df[!,:temp] .-273.15, df[!,:stability])
# ylabel!("Stability")
# p3 = plot(df[!,:temp] .-273.15, df[!,:biomass])
# ylabel!("Biomass")

# plot(p1, p2, p3, layout=(1,3))
