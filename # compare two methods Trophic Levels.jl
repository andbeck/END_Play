# compare two methods Trophic Levels

using EcologicalNetworks, EcologicalNetworksPlots
using Plots, Random, DataFrames, CSV
using LinearAlgebra

# Eva's Trophic Postion Code.

function normalize_matrix(A)
    A2 = transpose(A)
    colsum = sum(A2, dims = 1)
    colsum[colsum .== 0] .= 1
    normA = (A2'./vec(colsum))'
    return normA
end

function trophic_position(A)
    S = size(A,1)
    if S < 3
        return trophic_rank(A)
    else
        Mt = normalize_matrix(A)
        m = Int.(zeros(S,S))
        [m[i,i] = 1 for i in 1:S] #fill diag with 1
        detM = det(m .- Mt')
        if detM != 0
            tp = \(m .- Mt', repeat([1], S))
        else
            tmp = m
            for i in 1:9
                tmp = tmp * Mt' .+ m
            end
            tp = tmp * repeat([1], S)
        end
        return tp
    end
end



Random.seed!(123)
Fweb = N = nichemodel(50,0.05)
eltype(N)
Fweb_bin=adjacency(Fweb)

# first by trophic level with random postion within
I = initial(FoodwebInitialLayout, Fweb)
plot(I, Fweb)
scatter!(I, Fweb)

## Attempt to explore TL calcs in EN

a=trophic_level(Fweb)
b=distance_to_producer(Fweb, f = maximum)
c=trophic_position(Fweb_bin)

x=repeat(["S"],20)
y=collect(1:20)
convert(y, string)

x.*y

aa=DataFrame(a)
bb=DataFrame(b)
cc = DataFrame(TL = c)

CSV.write("TL.csv",aa)
CSV.write("d2p.csv",bb)
CSV.write("EvaTL.csv",cc)

tt=append!(aa, bb)
ttt = append!(tt, cc)
oo=append!(ttt, dd)

permutedims(oo)
