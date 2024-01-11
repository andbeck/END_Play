# compare two methods Trophic Levels

using EcologicalNetworks
using Plots, Random, DataFrames

Random.seed!(123)
Fweb = N = nichemodel(20,0.01)
eltype(N)

## Attempt to explore TL calcs in EN

a=trophic_level(Fweb)
b=distance_to_producer(Fweb, f = maximum)
c=distance_to_producer(Fweb, f = Statistics.mean)

aa=DataFrame(a)
bb=DataFrame(b)
cc=DataFrame(c)

tt=append!(aa, bb)
append!(tt, cc)


# steps in d2p function
Y = nodiagonal(Fweb)
# oops
paths = dijkstra(Y)
# works
paths = dijkstra(Fweb)

#
consumers = collect(keys(filter(p -> iszero(p.second), degree(Y; dims=1))))
filter!(int -> int.to ∈ consumers, paths)
tl = Dict{eltype(species(Y)),Float64}()
