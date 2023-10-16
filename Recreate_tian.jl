## Recreate Tians paper using their body sizes

## Packages
# Install
Pkg.activate(".")
Pkg.add("Plots")
Pkg.add("Statistics")
Pkg.add("Random")
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("StatsPlots")
Pkg.add("DifferentialEquations")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("DelimitedFiles")
Pkg.add("EcologicalNetworks")
Pkg.add("EcologicalNetworksPlots")

# Packages
using DataFrames, Plots, Random, Distributions
using EcologicalNetworks, EcologicalNetworksPlots, EcologicalNetworksDynamics
using LinearAlgebra, EcologicalNetworksDynamics

# set the seed
Random.seed!(12325)

function cascade_model_alina(C; mprod = [2, 5], minvert=[4, 6, 7, 8])
    nprod = length(mprod) #find the number of producers
    ninvert = length(minvert) #find the number of invertibrates
    phyto_names = "p".* [string(i) for i in 1:nprod] #give the phytoplankton names by looping over the number of producers and turn that number into a string (word) and then at P onto it by broadcasting p* over all the values from the loop 
    invert_names = "zoo".* [string(i) for i in 1:ninvert] #give the zooplankton names using the same method as above
    S = nprod + ninvert #calculate the species richness 
    idx_sort_mass_total = sortperm([mprod; minvert]) #sort the producer and invertibrate masses in order and gives a vector with the id (the position of the species in the orrigional vector), the ; in a square bracket means that it concatenates (joins) the two things together

    # Build the matrix of possible interaction
    M = ones(Bool,S, S) #make a matrix of 1s the size of the number of species, bool means that the 1s == true and 0s == false which is usefull later when some parts of the matrix are turned to 0 as this says if an interaction is there or not

    # Select where the link are authorized
    ## Select possible links: invertebrate - producers
    lower_triangle = tril(M, -1) #selects the lower triangle and so sets the upper triangle to 0, the -1 means it starts this from 1 row from the top so that the diagnonal is not selected
    ## Remove row of producers
    prod_idx = findall(idx_sort_mass_total .∈ [1:nprod]) #make a list of producer ids 
    invert_idx = findall(idx_sort_mass_total .∈ [nprod+1:S])
    idx = [prod_idx; invert_idx]

    possible_links = lower_triangle[setdiff(1:S, prod_idx),:] # : = all columns

    # Number Link in the foodweb
    L = C * (S-1) * S # actual number of links
    num_possible_links = sum(possible_links) #number of all possible links


    # Probability equal to the density of link:
    p = L / num_possible_links # the probability of a link forming is == to number of actual links/number of possible links

    if L > num_possible_links
        error("L > num of possible links, please decrease C!")
    end

    # Generate random uniform number
    sampling_link = rand(num_possible_links)
    # Set Links
    ## No interaction if random number sup or equal to p
    sampling_link[findall(>=(p), sampling_link)] .= 0
    ## Interaction for the remaining, i.e. the non zeros:
    sampling_link[findall(>(0), sampling_link)] .= 1

    possible_links[possible_links .== true] = sampling_link

    A = []
    n = 1
    for i in 1:S
        if i ∈ prod_idx
            push!(A, repeat([0], S))
        else
            push!(A, possible_links[n, :])
            n += 1
        end
    end

    A = reduce(hcat, A)'
    # Reorder columns as the provived vector of producer masses and zoomasses [mprod; minvert]
    A_resorted = A[idx, idx]
    fw = FoodWeb(A_resorted, M = [mprod; minvert], species = [phyto_names; invert_names],
                 metabolic_class = [repeat(["producer"], nprod); repeat(["invertebrate"], ninvert)]
    )
    fw
end

fw = cascade_model_alina(.2; mprod = [1, 10], minvert = [5, 8, 15])
# Check the final connectance
sum(fw.A) / ((richness(fw) - 1) * richness(fw)) ## all good