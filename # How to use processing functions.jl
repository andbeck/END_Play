# How to use processing functions

using Plots
using Random, DataFrames, StatsPlots
using EcologicalNetworks, EcologicalNetworksPlots
using EcologicalNetworksDynamics

Random.seed!(123)
S = 10
A = FoodWeb(nichemodel, S; C = 0.1)
params = ModelParameters(A)
bo = rand(S)
out = simulate(params, bo)

plot(out)


# biomass
biomass(out)
biomass(out).total

# diversity
EcologicalNetworksDynamics.richness(out)
species_persistence(out)
evenness(out)
shannon_diversity(out)

# stability
coefficient_of_variation(out)
community_cv(out)
species_cv(out)
