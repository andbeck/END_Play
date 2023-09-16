using EcologicalNetworksDynamics
using Random, DataFrames, Plots, CSV

######################
# vectors of variables
######################

#T_range = 0:2:40 # temperatures
T_range = 20:1:40
K_range = 1.:1:20 # K valus

# convert T to kelvin for BA
T_values = 273.15 .+ collect(T_range) # temperature 1-40C

# organis some K thing.s
K_int_values = collect(K_range) # intercept of the carrying capacity (eutrophication)
<<<<<<< Updated upstream
#n_T, n_K = length(T_range), length(K_range)
=======
>>>>>>> Stashed changes

# basal species starting biomass
m0 = 0.01


######################
# food webs from niche model
######################

# make 10 food web using niche model and Z value 100
Random.seed!(22)
FWs = []
for _ in 1:n_rep
    fw = FoodWeb(nichemodel, 30; C=0.1, Z=100)
    fw.M *= m0
    push!(FWs, fw)
end

# n vals for counting
n_T, n_K = length(T_range), length(K_range)
n_rep = 10 # reps
n_web = length(FWs)

# dataframe to store results
<<<<<<< Updated upstream
#n_lines = n_T * n_K * n_rep
df = DataFrame(step=[], temp=[], richness=[], stability=[], biomass=[])

# inital web
#foodweb = FWs[1]

    # declarations

simTemp = function (foodweb, tempRange)
    


    # initial paramters
    p = ModelParameters(foodweb,
        functional_response=ClassicResponse(foodweb, h=1.2),
        biorates=BioRates(foodweb; d=0))

=======
n_lines = n_web * n_T * n_K * n_rep
df = DataFrame(fw = [], step = [], temp=[], richness = [], stability = [], biomass = [])

for f in 1:n_web

    # start model
    # inital web
    foodweb = FWs[f]

    # initial paramters
    p = ModelParameters(foodweb,
                        functional_response=ClassicResponse(foodweb, h=1.2),
                        biorates=BioRates(foodweb; d=0))

>>>>>>> Stashed changes
    # set initial biomasses
    B0 = zeros(richness(foodweb))
    K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])
    B0[producers(foodweb)] .= K_prod
<<<<<<< Updated upstream
    B0[1:richness(foodweb).∉[producers(foodweb)]] .= K_prod / 8

    # run model across range of T's
    for i in 1:size(tempRange, 1)

        # make p = p_next for looping over T
        if i == 1
            p_next = p
            B0 = B0
            fw_next = foodweb
        else
            #p_next = p_next;
            B0 = biomass_next
            fw_next = A_next
        end

        TT = tempRange[i]
        KK = 1
        K_int = K_int_values[1]
        K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])

        # set temperature
        set_temperature!(p_next, TT, ExponentialBA(K=exp_ba_carrying_capacity(aₚ=K_int)))

        # simulate biomass dynamics for 10 years
        out = simulate(p_next, B0, tmax=3153600000,
            callback=EcologicalNetworksDynamics.ExtinctionCallback(1e-12, p_next, false),
            adaptive=true,
            dt=24 * 60 * 60,
            saveat=24 * 60 * 60,
        )

        # collect deets and write to df
        rr = richness(out)
        ss = community_cv(out)
        bb = biomass(out).total
        push!(df, [i, TT, rr, ss, bb])
        #df

        # identify extinctions and make mask
        # figure out how to deal with scenario with NO extinctions when who_extinct is empty
        who_extinct = keys(get_extinct_species(out))

        # a list of 1:n species with 0's in place of extinctions
        species = 1:size(fw_next.A, 1)
        #extant = setdiff(species, who_extinct)
        if !isempty(who_extinct)
            species_idx = setdiff(species, who_extinct) # or findall(x->x ∉ who_extinct, species)
        else
            species_idx = 1:size(fw_next.A, 1)
        end

        # define biomass vector for next iteration       
        biomass_next = biomass(out, last=1).species[species_idx]

        # define the food web for the next iteration
        # subset the matrix
        A_next = fw_next.A[species_idx, species_idx]
        # S-size(A_next)[1] == length(extinct_1) # a test
        A_next = FoodWeb(A_next)
        
        # define the bodymass vector for the next iteration
        A_next.M = fw_next.M[species_idx]
        #A_next

        # reset the params and bodymass vector with subsetted bodymass
        p_next = ModelParameters(A_next, functional_response=ClassicResponse(A_next, h=1.2), biorates=BioRates(A_next; d=0))

    end

    return df
end

outSim = simTemp(FWs[1], T_range)

p1 = plot(df[!, :temp] .- 273.15, df[!, :richness])
p2 = plot(df[!, :temp] .- 273.15, df[!, :stability])
p3 = plot(df[!, :temp] .- 273.15, df[!, :biomass])
=======
    B0[1:richness(foodweb) .∉ [producers(foodweb)]] .= K_prod / 8

    # run model across range of T's
    for i in 1:size(T_range,1)

        #################################################
        ## Organise things for temperature depdendence ##
        #################################################

        # define T for time step i
        TT = T_values[i]

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

# Visualisation
df
CSV.write("tempRun.csv", df)

>>>>>>> Stashed changes

# p1 = plot(df[!,:temp] .-273.15, df[!,:richness])
# ylabel!("Richness")
# p2 = plot(df[!,:temp] .-273.15, df[!,:stability])
# ylabel!("Stability")
# p3 = plot(df[!,:temp] .-273.15, df[!,:biomass])
# ylabel!("Biomass")

# plot(p1, p2, p3, layout=(1,3))
