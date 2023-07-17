# 14.07.2023
# extinction of apex predator

# Packages
using Random, Plots, Distributions, DataFrames, StatsPlots
using EcologicalNetworksDynamics

# Initial Food Web
fw = FoodWeb([9 => [1,2],
10 => [1,2,3,4,5,6],
11 => [5,6],
12 => [5,6],
13 => [7,8],
14 => [9],
15 => [9,10],
16 => [12,13],
17 => [10],
18 => [10,11,13]], Z = 10)

# See the matrix
transpose(fw.A)
# Setting biomass of 18 species (random)
B0 = rand(18)
# Setting model parameters
params = ModelParameters(fw)




# Burn in
out = simulate(params, B0)
plot(out) # good = equilibrium
# Note the extinctions that happened
get_extinct_species(out)
# Community at the end of burn in
coefficient_of_variation(out)
richness(out)
# Look at web at end of burn in
fw_2 = trophic_structure(out).alive_A
trophic_classes(fw_2)



# Extra Information
trophic_classes(fw)
trophic_levels(fw)
trophic_structure(out)

# overall biomass by community and species
bb = biomass(out)
# can calcuate average biomass of producers etc using information
# from trophic_classes
bb.total
# just producers
bb.species[1:8]
# or using the indices from trophic_classes!!
bb.species[trophic_classes(fw).producers]






# Calculate biomass
biomass(out)
# Biomass in %
species_persistence(out)
# Make extinction happen
# identify apex predator at the end of the burn in
apex = findmax(trophic_levels(fw_2))
apex[2] # this is the species that is at the top










# Deepcopy FW
fw_3 = deepcopy(fw_2)

# Create food web with deliberate extinction

# Step 1: create list of numbers but missing the index/id of the apex predator
# .!=apex is 'broadcast to all the number the desire to leave out the apex identity (no. 14)
mask_topPred = 1:size(fw_3,1) .!=apex[14]



# Step 2: use this `mask`` to subset the matrix
# we use the [mask, mask] syntax for both rows and columns to be left out.
fw_3 = fw_3[mask_topPred, mask_topPred]
# check
size(fw_2)
size(fw_3)

# Step 4: create a set of indices where biomasses are non-zero (e.g. the not extinct species)
mask_biom = biomass(out, last = 1).species .!=0

# Step 5: use this mask to isolate the final biomasses of the non-extinct
# species in the burn in.
B0_end = biomass(out, last = 1).species[mask_biom]

# Step 6: remove the apex predator from this list
# you are making it extinct here!
B0_restart = B0_end[mask_topPred]

### RESTARTING THE MODELLING ####
# Step 1: create the new food web
fw_restart = FoodWeb(fw_3)

# Step 2: re-initialise the parameters
params_restart = ModelParameters(fw_restart)

# Step 3: re-simulate.
out_end = simulate(params_restart, B0_restart, tmax = 1000)

#### PLOTTING THE TRAJECTORIES ####
# rather complicated, but it works and will explain later!
plot(out, c = :grey50, leg = false, xlims=(0, out_init.t[end]+out_end.t[end]))
# add the second sim... note different format to allow for 
# re-defining start time to end of initial series.
plot!(out_end.t .+ out_init.t[end], out_end', c = :red)
