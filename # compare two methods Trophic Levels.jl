# compare two methods Trophic Levels

using EcologicalNetworks
using Plots, Random, DataFrames

Random.seed!(123)
Fweb = N = nichemodel(20,0.1)
eltype(N)

# first by trophic level with random postion within
I = initial(FoodwebInitialLayout, Fweb)
plot(I, Fweb)
scatter!(I, Fweb)

## Attempt to explore TL calcs in EN

a=trophic_level(Fweb)
b=distance_to_producer(Fweb, f = maximum)
c=distance_to_producer(Fweb, f = Statistics.mean)
d = distance_to_producer(Fweb, f = minimum)

aa=DataFrame(a)
bb=DataFrame(b)
cc=DataFrame(c)
dd=DataFrame(d)

tt=append!(aa, bb)
ttt = append!(tt, cc)
oo=append!(ttt, dd)

permutedims(oo)
