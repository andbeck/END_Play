using EcologicalNetworksDynamics, Plots, Random

Random.seed!(123)

#### The INITIALISATION ####

# first create a food web of 20 species
# 5% connectance and PPMR of 10x
foodweb = FoodWeb(nichemodel, 20, C = 0.05, Z = 10)

# set the biomasses of the 20 species
B0 = rand(20)

# set up the model parameters
# could do more here with the functional response etc
params = ModelParameters(foodweb)

#### The BURN IN PHASE ####

# create BURN IN PHASE
out_init = simulate(params, B0, verbose = false)

# Could do some processing to find out what the community looks like and behaved like

# Note the extinctions that happend
get_extinct_species(out_init)

cv_burn = coefficient_of_variation(out_init)
rich_burn = richness(out_init)
ts_burn = trophic_structure(out_init).alive_A
tc_burn = trophic_classes(ts_burn)

#############################################
### identify extinctions in prep for i+1 ####
#############################################

    # who is extinct
    who_extinct = keys(get_extinct_species(out_init))

    # a list of 1:n species in current network
    species = 1:size(foodweb.A, 1)

    # mask extinct species or leave as is if who_extinct is empty
    if !isempty(who_extinct)
        species_idx = setdiff(species, who_extinct) # or findall(x->x ∉ who_extinct, species)
    else
        species_idx = 1:size(foodweb.A, 1)
    end

    species_idx
######################################
## Update foodweb, B0 and p for t+1 ##
######################################

    # subset the matrix
    foodweb = FoodWeb(foodweb.A[species_idx, species_idx])

    # subset and collect the biomass (last value approach vs. mean?)
    # whether the mean or the last value is taken doesn't seem to matter
    B0 = biomass(out_init, last = 1).species[species_idx]

    ## COULD DEFINE NEW VALUES OF B0 here... e.g. harvesting of a species. ###

    # reset the params and bodymass vector with subsetted bodymass
    p = ModelParameters(foodweb)

    out_next = simulate(p, B0, tmax = 2000, verbose = false)

# plotting
# note we have out_init and out_next

# first plot the burnin (out_init), but set the x limits to max of 
# the first + second simulation (out_next)
plot(out_init, c = :grey50, leg = false, 
    xlims=(0, out_init.t[end]+out_next.t[end]))

# now add the second sim, starting at the end of the burn-in
# note that it is an apostrophe after out_next, it is not a typo.

plot!(out_next.t .+ out_init.t[end], 
    out_next') n
    