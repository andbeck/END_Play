using EcologicalNetworksDynamics, Plots, Random

Random.seed!(123)

#### The INITIALISATION ####

# first create a food web of 20 species
# 5% connectance and PPMR of 10x
fw = FoodWeb(nichemodel, 20, C = 0.05, Z = 10)

# set the biomasses of the 20 species
B0 = rand(20)

# set up the model parameters
# could do more here with the functional response etc
params = ModelParameters(fw)

#### The BURN IN PHASE ####

# create BURN IN PHASE
out_init = simulate(params, B0, tmax = 2000, verbose = false)

# Note the extinctions that happend
get_extinct_species(out_init)

# Could do some processing to find out what the community looks like and behaved like
coefficient_of_variation(out_init)
richness(out_init)

# Look at the food web at end of burn in
fw_2 = trophic_structure(out_init).alive_A
trophic_classes(fw_2)

# identify the apex predator at the end of the burn in.
# use findmax and our trophic_levels function
apex = findmax(trophic_levels(fw_2))
apex[2] # this is the species that is at the top

# THIS IS V. IMPORTANT
# if you don't do this, everything from here will over-write fw_2, which we want to keep.
fw_3 = deepcopy(fw_2)

###  CREATING THE FOOD WEB WITH A DELIBERATE EXTINCTION ###

# Step 1: create list of numbers but missing the index/id of the apex predator
# 1:size(fw_3, 1) is a sequence of numbers == S species (not extinct)
# .!=apex is 'broadcast to all the number the desire to leave out the apex identity (no. 11)
mask_topPred = 1:size(fw_3,1) .!=apex[2]

# Step 2: use this `mask`` to subset the matrix
# we use the [mask, mask] syntax for both rows and columns to be left out.
fw_3 = fw_3[mask_topPred, mask_topPred]

# Step 3: check it worked.  
size(fw_2)
size(fw_3)

# Step 4: create a set of indices where biomasses are non-zero (e.g. the not extinct species)
mask_biom = biomass(out_init, last = 1).species .!=0

# Step 5: use this mask to isolate the final biomasses of the non-extinct
# species in the burn in.
B0_end = biomass(out_init, last = 1).species[mask_biom]

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
plot(out_init, c = :grey50, leg = false, xlims=(0, out_init.t[end]+out_end.t[end]))
# add the second sim... note different format to allow for 
# re-defining start time to end of initial series.
plot!(out_end.t .+ out_init.t[end], out_end', c = :red)