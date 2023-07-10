# make invasions happen

using EcologicalNetworksDynamics, Plots, Random
using StatsBase

Random.seed!(12234)

#### The INITIALISATION ####

S = 20

# first create a food web of 20 species
# 5% connectance and PPMR of 10x
fw = FoodWeb(nichemodel, S, C = 0.05, Z = 10)

##### second - identify an invader ######
# we will use a species already in the community as the invader

# Option 1: top predator
# identify the apex predator in the food web
# use findmax and our trophic_levels function
apex = findmax(trophic_levels(fw))
intro_apex = apex[2]

# Option 2: an intermediate consumer
int = trophic_classes(fw).intermediate_consumers
intro_intermediate = sample(int)

#### Create FW without Invader
# we do this by setting it's biomass to 0.
# this means that the species has a placeholder in the foodweb,
# but won't influence the burn-in dynamics

# set the biomasses of the 20 species
B0 = rand(20)
B0[intro_apex] = 0 # set the 19th element to 0
B0 # check it

# Build the Burn In

# could do more here with the functional response etc
params = ModelParameters(fw)
out_init = simulate(params, B0)

# it is possible that lots of species go extinct
# we want to remove networks where <80 of the species persist
# persistence = species_persistence(out_init, last = 1)

# Collect Data on Burnin
init_persistence = species_persistence(out_init, last = 1)
init_biomass = biomass(out_init, last = 2)
init_stability = community_cv(out_init)

### MANAGEMENT OF INVADER and EXTINCTIONS
# the invasive species we set at the start will be in the list

  # Collect the Deterministic Extinctions from Burnin ----

  # Step 1: get all the extinctions
  # note that the one we set to 0 (19) is in the list
  # note that the 0.0 in this this object is the TIME
  # that the species went extinct.
  extinctions = get_extinct_species(out_init)

  # remove the introduced ID from the extinction tracker
  # this is a tricky bit of code, but basiucally it's saying
  # keep only the non-zero things in the second column of the exctinctions thingy
  filter!(p -> !=(0.0, p.second), extinctions)

  # check it
  extinctions

  ## NOW we set introduced biomass of the invader to be > 0 ----

  # Step 1: get the biomasses at the end of the first simulation above
  B0_Inv = out_init.u[end,:][1]

# step 2: set the introduced to a positive biomass
# we'll use the mean biomass of the original species

  B0_Inv[intro_apex] = mean(B0)

# now simulate again starting with the B0_Inv numbers
# note that the invader is at a positive biomass,
# but than any extinct species are at 0
  out_invade = simulate(params, B0_Inv)



#### PLOTTING THE TRAJECTORIES ####

lty = fill(:solid, 1, S) #original species: solid line
lty[intro_apex] = :dash #invader: dash line
col = fill(:red, 1, S)
col[intro_apex] = :green

# rather complicated, but it works and will explain later!
plot(out_init, c = :grey50, leg = false, xlims=(0, out_init.t[end]+out_invade.t[end]))
# add the second sim... note different format to allow for 
# re-defining start time to end of initial series.
plot!(out_invade.t .+ out_init.t[end], out_invade', c = col, linestyle = lty)
scatter!([out_init.t[end]], [B0_Inv[intro_apex]], mc = :black, ms = 5)
