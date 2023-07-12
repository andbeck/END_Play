# make invasions happen with specific targets

using EcologicalNetworksDynamics, Plots, Random
using StatsBase

Random.seed!(123)

# make a food web
fw = FoodWeb([8=>[7,6,5], 7=>[4,3], 6=>[3,1], 5=>[3,2]])
fw.A

trophic_classes(fw)

# we want to insert a an intermediate consumer that preys on
# 2/3 of the resources but is NOT eaten by a predator

resources = trophic_classes(fw).producers
whichresources = sample(resources, 3, replace = false)

fw2 = deepcopy(fw)

# the predation matrix
A = Matrix(fw2.A)

# the new species (eats which resources)
intSpRow = zeros(Int8, 8)
intSpRow[whichresources] .= 1

# the new species as a resource (eaten by none)
intSpCol= zeros(Int8, 9)

AA = [A[]]


# build the matrix
A2 = vcat(A, intSpRow')
A2 = hcat(A2, intSpCol)

# make it a food web
fw3 = FoodWeb(A2)
trophic_classes(fw3)

(A2[1:7,:], A2[9,:], A2[8,:])