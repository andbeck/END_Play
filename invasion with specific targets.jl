# make invasions happen with specific targets

using EcologicalNetworksDynamics, Plots, Random
using StatsBase

Random.seed!(123)

# define 'no beaver' food web
initLinks = [8=>[7,6,5], 7=>[4,3], 6=>[3,1], 5=>[3,2]]
initFW = FoodWeb(initLinks)

# find all the trees
resources = trophic_classes(initFW).producers
# select a set of the trees
whichresources = sample(resources, 3, replace = false)

# Define what the beaver eats
# and that it is not eaten
newLinks = [9=>whichresources]

# combine these new links with the initial links
invadedLinks = vcat(initLinks, newLinks)

# create the invaded by beaver food web
invadedFW = FoodWeb(invadedLinks)

# Look at the matrices if you want
initFW.A
invadedFW.A

# Run the burn-in
initParams = ModelParameters(initFW)
initB0 = rand(8)
initSol = simulate(initParams, initB0)

# set up the invasion biomasses.
# end of sims biomass for original 8,
# and mean of these for beaver.
invadedB0 = vcat(biomass(initSol, last = 1).species,
            mean(biomass(initSol, last = 1).species))

#
invadedParams = ModelParameters(invadedFW)
invadedSol = simulate(invadedParams, invadedB0)

plot(initSol, c = :grey50, leg = false,
        xlims=(0, initSol.t[end]+invadedSol.t[end]))

col = fill(:red, 1, 9)
col[9] = :green

plot!(invadedSol.t .+ initSol.t[end], invadedSol', c = col)
scatter!([initSol.t[end]], [invadedB0[9]])
