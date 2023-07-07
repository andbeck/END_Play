using EcologicalNetworksDynamics, Plots, Random

Random.seed!(123)

fw = FoodWeb(nichemodel, 20, C = 0.05, Z = 10)
B0 = rand(20)

params = ModelParameters(fw)

out_init = simulate(params, B0, tmax = 2000, verbose = false)
get_extinct_species(out_init)

# find matrix at end of burn - in
fw_2 = trophic_structure(out_init).alive_A
trophic_classes(fw_2)

# identify apex
apex = findmax(trophic_levels(fw_2))
apex[2] # this is the species that is at the top

# THIS IS V. IMPORTANT
fw_3 = deepcopy(fw_2)

# Now create list of numbers missing the apex
mask_topPred = 1:size(fw_3,1) .!=apex[2]

# use this mask to subset the matrix
fw_3 = fw_3[mask_topPred,mask_topPred]
size(fw_3)

mask_biom = biomass(out_init, last = 1).species .!=0
B0_end = biomass(out_init, last = 1).species[mask_biom]
B0_restart = B0_end[mask_topPred]

# restart with apex predator missing
fw_restart = FoodWeb(fw_3)
params_restart = ModelParameters(fw_restart)
out_end = simulate(params_restart, B0_restart)
