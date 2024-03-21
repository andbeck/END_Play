using EcologicalNetworksDynamics, Plots, Random
using DataFrames

#### The INITIALISATION ####

# # first create a food web of 20 species
# # 5% connectance and PPMR of 10x
# foodweb = FoodWeb(nichemodel, 20, C = 0.05, Z = 10)

foodweb = FoodWeb(
[6 => [1,2,3,4,5],
7 => [1,2,3,4,5,6],
8 => [1,2,3,4,5,6],
9 => [1,2,3,4,5],
10 => [1,2,3,4,5,6],
11 => [1,2,3,4,5,6],
12 => [1,2,3,4,5,6],
13 => [7,8,9,10,11,12],
14 => [6,7,8,9,10,11,12],
15 => [7,8,9,10,11,12,13,14]], Z = 10)

#randomise starter
Random.seed!(123)

# INTIAL SETUP
# set the biomasses of the 20 species
# parrot is 12?
B0 = rand(15)

# set up the model parameters
# could do more here with the functional response etc
params = ModelParameters(foodweb)


## Collection Zone - set up data frame to collect richness, biomass and stability
df = DataFrame(richness = [], biomass = [], stability = [], Sp12_Biomass = [])

## Looping

for i in 1:10

    # simulate -> i = 1 is burn-in
    out = simulate(params, B0, verbose = false)

    # collect data for each cycle
    push!(df, [richness(out), biomass(out).total, coefficient_of_variation(out).community,
        B0[12]])

    #############################################
    ### identify extinctions in prep for i+1 ####
    #############################################

    # who is extinct
    who_extinct = keys(get_extinct_species(out))

    # a list of 1:n species in current network
    species = 1:size(foodweb.A, 1)

    # mask extinct species or leave as is if who_extinct is empty
    if !isempty(who_extinct)
        species_idx = setdiff(species, who_extinct) # or findall(x->x ∉ who_extinct, species)
    else
        species_idx = 1:size(foodweb.A, 1)
    end


    @info "Iteration number: $i"
    @info "Extinction: $who_extinct"
    B12 = B0[12]
    @info "Biomass Sp 12: $B12"

    ######################################
    ## Update foodweb, B0 and p for t+1 ##
    ######################################

    # subset the matrix
    foodweb = FoodWeb(foodweb.A[species_idx, species_idx])

    # subset and collect the biomass (last value approach vs. mean?)
    # whether the mean or the last value is taken doesn't seem to matter
    B0 = biomass(out, last = 1).species[species_idx]

    #### reduce 12th species biomass by 10% ####
    #### doesn't actually work, as species 12 is not the same after extinctions
    #### how would you manage this?

    B0[12] = B0[12]*0.90

    # reset the params and bodymass vector with subsetted bodymass
    params = ModelParameters(foodweb)
end

# see the data
df
df[!,:Sp12_Biomass]

# plotting
step = 1:10
a = plot(step, df[!,:Sp12_Biomass])
title!("Species12")
b = plot(step, df[!,:biomass])
title!("Biomass")
c = plot(step, df[!,:stability])
title!("Stability")
d = plot(step, df[!,:richness])
title!("Richness")

plot(a,b,c,d, layout=(2,2))
