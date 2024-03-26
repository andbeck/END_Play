using EcologicalNetworksDynamics, Plots, Random
using DataFrames

#### The INITIALISATION ####

# # first create a food web of 20 species
# # 5% connectance and PPMR of 10x
# foodweb = FoodWeb(nichemodel, 20, C = 0.05, Z = 10)

# The parrot community food web
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
Random.seed!(1234)

# INTIAL SETUP
# set the biomasses of the 20 species
# parrot is 8?

B0 = rand(15)

# set up the model parameters
# can add Mass (body size) details here
# could do more here with the functional response etc
# (e.g. default is a type 3 functional response, which is very stable)
# (what happens when it's a type 2 and less stable)

params = ModelParameters(foodweb)

## Collection Zone - set up data frame to collect richness, biomass and stability
df = DataFrame(richness = [], biomass = [], stability = [], Sp8_Biomass = [])

# collect all species biomasses at end of each simulation
all_biomass = Any[]

## Looping

# 10 harvesting events
for i in 1:10

    # just some information during the simulations
    B8_start = B0[8]
    @info "Iteration number: $i"
    @info "Biomass Sp 8 init: $B8_start"

    # simulate
    out = simulate(params, B0, verbose = false)

    # collect data for each cycle
    push!(df, [richness(out), biomass(out).total, coefficient_of_variation(out).community,
        B0[8]])
    push!(all_biomass, B0)

    ######################################
    ## NEW: Update Food Web with Zero's for extinct species
    ## Retains structure and indexing for plotting
    ######################################

    # grab the trophic structure
    troph_struc = trophic_structure(out, last = 1)

    # grab the alive species
    alive_species = troph_struc.alive_species

    # pre-allocate new biomass with all zeros
    B0 = zeros(size(foodweb.A)[1])

    # replace with biomass of alive species at end of sim
    # any that are extinct are assigned a 0.
    B0[alive_species] = biomass(out, last = 1).species[alive_species]

    # just some information during the simulations
    B8_end = B0[8]
    @info "Biomass Sp 8 end: $B8_end"


    # harvest parrot at 10%
    B0[8] = B0[8]*0.9

    # just some information during the simulations
    B8_harvest = B0[8]
    @info "Biomass Sp 8 harvest: $B8_harvest"

    # reset the params and bodymass vector with subsetted bodymass
    # note no change to the food web structure as we simply mark
    # extinct ones as 0 biomass
    params = ModelParameters(foodweb)

end

# see the data frame
df

# see parrot trajectory
df[!,:Sp8_Biomass]

## convert all the B0 vectors (all species) to a data frame
mat = zeros(length(all_biomass), length(all_biomass[1])) # create a matrix the right size

# stick the data in the right place
for i in 1:length(all_biomass)
    mat[i,:] = all_biomass[i]
end

# see the matrix
# each column is a species biomass through time.
mat

# OR classy!
#DataFrame(mapreduce(permutedims, vcat, all_biomass[:,:]), :auto)

# plotting
step = 1:10
a = plot(step, df[!,:Sp8_Biomass])
title!("Species8")
b = plot(step, df[!,:biomass])
title!("Biomass")
c = plot(step, df[!,:stability])
title!("Stability")
d = plot(step, df[!,:richness])
title!("Richness")
e = plot(step, mat, legend = false)
title!("All Species ~ Time")

plot(a,b,c,d,e, layout=(2,3))
