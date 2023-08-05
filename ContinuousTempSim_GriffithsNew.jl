using EcologicalNetworksDynamics
using Random, DataFrames

# vectors of variables
T_range = 0:2:40
K_range = 1.:1:20
n_rep = 10
T_values = 273.15 .+ collect(T_range) # temperature 1-40C
K_int_values = collect(K_range) # intercept of the carrying capacity (eutrophication)
n_T, n_K = length(T_range), length(K_range)

# basal species starting biomass
m0 = 0.01

# make 10 food web using niche model and Z value 100
Random.seed!(22)
FWs = []
for _ in 1:n_rep
    fw = FoodWeb(nichemodel, 30; C=0.1, Z=100)
    fw.M *= m0
    push!(FWs, fw)
end


# dataframe to store results
n_lines = n_T * n_K * n_rep
df = DataFrame(step = [], temp=[], richness = [])

# start model
foodweb = FWs[1]
p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2), biorates=BioRates(foodweb; d=0))
# set initial biomasses
B0 = zeros(richness(foodweb))
K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])
B0[producers(foodweb)] .= K_prod
B0[1:richness(foodweb) .∉ [producers(foodweb)]] .= K_prod / 8


#for i in 1:2
    TT = T_values[i]
    KK = 1
    K_int = K_int_values[1]
    K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])

    # set temperature
    set_temperature!(p, TT, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K_int)))

    # simulate biomass dynamics for 10 years
    out = simulate(p, B0 ,tmax = 3153600000,
    callback = EcologicalNetworksDynamics.ExtinctionCallback(1e-12, p, false),
    adaptive = true,
    dt = 24*60*60,
    saveat = 24*60*60,
    )

    # collect deets and write to df
    rr = richness(out)
    push!(df, [i, TT, rr])
    df
    # update Biomass
    B0 = biomass(out, last = 1).species

    # update foodweb
    extinctions = get_extinct_species(out)
    whichDel = keys(extinctions)
    fwA2 = trophic_structure(out).alive_A

    species_keep = 1:length(B0) .!= whichDel

    # Update the food web.... HOW??!!
    foodweb.A = fwA2
    foodweb.M[species_keep]
    p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2), biorates=BioRates(foodweb; d=0))
end
