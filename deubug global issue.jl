using EcologicalNetworksDynamics
using Random, DataFrames, Plots

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
df = DataFrame(step = [], temp=[], richness = [], stability = [], biomass = [])

# start model
# inital web
foodweb = FWs[1]

# initial paramters
p = ModelParameters(foodweb,
                    functional_response=ClassicResponse(foodweb, h=1.2),
                    biorates=BioRates(foodweb; d=0))

# set initial biomasses
B0 = zeros(richness(foodweb))
K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])
B0[producers(foodweb)] .= K_prod
B0[1:richness(foodweb) .∉ [producers(foodweb)]] .= K_prod / 8

# run model across range of T's
for i in 1:size(T_range,1)

    # make p = p_next for looping over T
    if i == 1
        p_next = p;
        B0 = B0;
        fw_next = foodweb;
    else
        #p_next = p_next;
        B0 = biomass_next;
        fw_next = A_next;
    end

    TT = T_values[i]
    KK = 1
    K_int = K_int_values[1]
    K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])

    # set temperature
    set_temperature!(p_next, TT, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K_int)))

    # simulate biomass dynamics for 10 years
    out = simulate(p_next, B0 ,tmax = 3153600000,
    callback = EcologicalNetworksDynamics.ExtinctionCallback(1e-12, p_next, false),
    adaptive = true,
    dt = 24*60*60,
    saveat = 24*60*60,
    )

    # collect deets and write to df
    rr = richness(out)
    ss = community_cv(out)
    bb = biomass(out).total
    push!(df, [i, TT, rr, ss, bb])
    #df

    # identify extinctions and make mask
    # figure out how to deal with scenario with NO extinctions when who_extinct is empty
    who_extinct = keys(get_extinct_species(out))

    # a list of 1:n species with 0's in place of extinctions
    species = 1:size(fw_next.A, 1)
    #extant = setdiff(species, who_extinct)
    if !isempty(who_extinct)
        species_idx = setdiff(species, who_extinct) # or findall(x->x ∉ who_extinct, species)
    else
        species_idx = 1:size(fw_next.A, 1)
    end

    # subset the matrix
    A_next = fw_next.A[species_idx, species_idx]
    # S-size(A_next)[1] == length(extinct_1) # a test
    A_next = FoodWeb(A_next)
    # subset the bodymasses
    A_next.M = fw_next.M[species_idx]
    #A_next

    # subset and collect the biomass (last value approach vs. mean?)
    # whether the mean or the last value is taken doesn't seem to matter
    biomass_next = biomass(out, last = 1).species[species_idx]

    # reset the params and bodymass vector with subsetted bodymass
    p_next = ModelParameters(A_next, functional_response=ClassicResponse(A_next, h=1.2), biorates=BioRates(A_next; d=0))
end

df

p1 = plot(df[!,:temp] .-273.15, df[!,:richness])
p2 = plot(df[!,:temp] .-273.15, df[!,:stability])
p3 = plot(df[!,:temp] .-273.15, df[!,:biomass])

plot(p1, p2, p3)
