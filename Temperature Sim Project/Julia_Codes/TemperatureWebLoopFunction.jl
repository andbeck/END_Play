##############################################
## MASTER FUNCTION TO DO THE CONTINGENCY WORK
##############################################
simTemp = function(webs, tempSeq)

    # set collection
    df = DataFrame(fw = [], step = [], temp=[], richness = [], stability = [], biomass = [])

    # over replicate (25) webs
    for f in 1:length(webs)

        # start model
        # inital web
        foodweb = webs[f]

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
        # from specific tempSeq

        for i in 1:size(tempSeq,1)

            #################################################
            ## Organise things for temperature depdendence ##
            #################################################

            # define T for time step i
            TT = tempSeq[i]+273.15 # adjust to Kelvin

            # set K vals
            KK = 1
            K_int = K_int_values[1]
            K_prod = unique(p.producer_growth.K[.!isnothing.(p.producer_growth.K)])

            # allocate temperature dependence to rates
            set_temperature!(p, TT, ExponentialBA(K = exp_ba_carrying_capacity(aₚ = K_int)))

            ## simulate biomass dynamics for 10 years ##
            out = simulate(p, B0 ,tmax = 3153600000,
                callback = EcologicalNetworksDynamics.ExtinctionCallback(1e-12, p, false),
                adaptive = true,
                dt = 24*60*60,
                saveat = 24*60*60,
            )

            ## collect deets and write to df ##
            rr = richness(out)
            ss = community_cv(out)
            bb = biomass(out).total
            push!(df, [f, i, TT, rr, ss, bb])

            #############################################
            ### identify extinctions in prep for i+1 ####
            #############################################

            # who is extinct
            who_extinct = keys(get_extinct_species(out))

            # a list of 1:n species in current network
            species = 1:size(foodweb.A, 1)
            #extant = setdiff(species, who_extinct)

            # mask extinct species or leave as is if who_extinct is empty
            if !isempty(who_extinct)
                species_idx = setdiff(species, who_extinct) # or findall(x->x ∉ who_extinct, species)
            else
                species_idx = 1:size(foodweb.A, 1)
            end

            ######################################
            ## Update foodweb, B0 and p for t+1 ##
            ######################################

            # subset the matrix
            foodweb = FoodWeb(foodweb.A[species_idx, species_idx])

            # subset and collect the biomass (last value approach vs. mean?)
            # whether the mean or the last value is taken doesn't seem to matter
            B0 = biomass(out, last = 1).species[species_idx]

            # reset the params and bodymass vector with subsetted bodymass
            p = ModelParameters(foodweb, functional_response=ClassicResponse(foodweb, h=1.2), biorates=BioRates(foodweb; d=0))
        end
    end
    return(df)
end
