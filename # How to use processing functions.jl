## Packages -------------

using Plots
using Random, DataFrames, StatsPlots
using EcologicalNetworksDynamics

## Repeatability/Sharing ----
Random.seed!(123)

## CREATING THE COMMUNITY ----------------
# Species Richness of the Food Web
S =20

# Create a Food Web using the Niche Model
A = FoodWeb(nichemodel, S; C = 0.1)

# Create a Food Web manually
A2 = FoodWeb([3=>(2,1)])

## Set up for the modelling -----------------
# Define parameters based on body size
params = ModelParameters(A)

# Set initial biomasses of each species
bo = rand(S)

# create simulation
out = simulate(params, bo)

# visualise
plot(out)

## POST PROCESSING -----------

# biomass (species, community)
biomass(out)
biomass(out).total

# diversity (richness, persistence, metrics of diversity)
richness(out)
species_persistence(out)
evenness(out)
shannon_diversity(out)

# stability (variability [others possible])
coefficient_of_variation(out)
community_cv(out)
species_cv(out)
