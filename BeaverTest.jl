using EcologicalNetworksDynamics, Random, Plots, StatsPlots, DataFrames
using StatsBase, CSV

#### Create data collection DataFrames
    biomass_df = DataFrame(repeat = [], iteration = [], ecosystem = [], trophic_class = [], biomass = [])
    stability_df = DataFrame(repeat = [], iteration = [], stability = [])
    shannon_df = DataFrame(repeat = [], iteration = [], shannon = [])

    ## and a quick dataframe checking if species don't get picked

#### CREATE METACOMMUNITY

    # max attempts
    max_attempts_meta = 1
    attempts = 1
    while true
        ## Create terrestrial foodweb
            global S_t = 6
            global C_t = 0.49
            global fw_t = FoodWeb(nichemodel, S_t, C = C_t)

        ## Create aquatic foodweb
            global S_a = 6
            global C_a = 0.49
            global fw_a = FoodWeb(nichemodel, S_a, C = C_a)

        ## Combine terrestrial and aquatic foodweb
            fw_prelinks = FoodWeb([fw_t.A zeros(S_a, S_t); zeros(S_t, S_a) fw_a.A])

        ## Create cross-ecosystem links

            # Define each trophic class per habitat

                #terrestrial
                global t_prod = trophic_classes(fw_prelinks).producers
                    filter!(x -> x <= S_t, t_prod)
                global t_int = trophic_classes(fw_prelinks).intermediate_consumers
                    filter!(x -> x <= S_t, t_int)
                global t_pred = trophic_classes(fw_prelinks).top_predators
                    filter!(x -> x <= S_t, t_pred)

                # aquatic
                global aq_prod = trophic_classes(fw_prelinks).producers
                    filter!(x -> x >= (S_t+1), aq_prod)
                global aq_int = trophic_classes(fw_prelinks).intermediate_consumers
                    filter!(x -> x >= (S_t+1), aq_int)
                global aq_pred = trophic_classes(fw_prelinks).top_predators
                    filter!(x -> x >= (S_t+1), aq_pred)

            # copy prelinks to post links
                fw_postlinks = deepcopy(fw_prelinks)

            # turn fw_postlinks into a matrix to add links
                global fw_postlinks_matrix = Matrix(fw_postlinks.A)

            # Create cross ecosystem links

                # First, terrestrial int and pred eating aquatic int and pred

                    # choose how many links to create
                        ta_links = round(S_t/5)

                    # Generate links from randomly sampled predators and prey
                        for _ in 1:ta_links
                            terrestrial = rand(vcat(t_int, t_pred))
                            aquatic = rand(vcat(aq_int, aq_pred))
                            
                            # Create links
                            fw_postlinks_matrix[terrestrial, aquatic] = 1
                        end

                # Then aquatic int and pred eating terrestrial int ONLY

                    # choose how many links to create
                        at_links = round(S_a/6)

                    # Generate links from randomly sampled predators and prey
                        for _ in 1:at_links
                            aquatic = rand(vcat(aq_int, aq_pred))
                            terrestrial = rand(t_int)
                            
                            # Create links
                            fw_postlinks_matrix[aquatic, terrestrial] = 1
                        end

            # Turn modified matrix back into foodweb object
                global metacom = FoodWeb(fw_postlinks_matrix, Z = 10) # Also give species bodymass ratio

        ## Create reference biomass values
            
            # generate parameters
            params = ModelParameters(metacom)
            
            # Set biomasses randomly between 0 and 1
            B_reference = rand(length(metacom.species))

            # run metacom until no species go extinct. Use these biomass values for the succession loop
            global sim_reference = simulate(params, B_reference, verbose = false)
            extinct_species = get_extinct_species(sim_reference)

            # ensure loop stops after max_attempts_meta
            attempts += 1
            if attempts > max_attempts_meta
                @warn "No suitable colonizer found after $max_attempts_meta attempts. Exiting the loop."
                break
            end


            if isempty(extinct_species) && all(!isempty(v) for v in [t_prod, t_int, t_pred, aq_prod, aq_int, aq_pred])
                break
            end
    end
    # check metacom
        metacom.A

    # Get equilibrium biomass values
        B_reference = biomass(sim_reference).species

    ## Create trophic class vectors for future indexing
        producers = trophic_classes(metacom).producers
        int_consumers = trophic_classes(metacom).intermediate_consumers
        predators = trophic_classes(metacom).top_predators
        metaspecies = vcat(producers, int_consumers, predators)

#### CREATE STARTING COMMUNITY

    # Potential issue: is there a way of ensuring the sampled t_int to be a beaver is a herbivore and/or has n links to producers?
    # give a tiny number of aquatic producers biomass?
    
    ## Define terrestrial and aquatic species, and a beaver
        terrestrial_species = vcat(t_prod, t_int, t_pred)
        aquatic_species = vcat(aq_prod, aq_int, aq_pred)
        beaver = sample(t_int)

    ## Set biomasses

        # Start with reference biomasses
        Bstart = deepcopy(B_reference)

        # Set aquatic and beaver species biomasses to 0
        Bstart[aquatic_species] .= 0
        Bstart[beaver] = 0

#### START REPEAT LOOP FOR STATISTICS
    for j in 1:1

        #### RUN STARTING COMMUNITY

            ## Generate parameters
                params = ModelParameters(metacom)

            ## Run starting community
                sim_start = simulate(params, Bstart)

            ## Collect data

                # Biomass
                    sim_start_bio_prod_t = sum(biomass(sim_start).species[t_prod])
                    sim_start_bio_int_t = sum(biomass(sim_start).species[t_int])
                    sim_start_bio_pred_t = sum(biomass(sim_start).species[t_pred])
                    sim_start_bio_prod_aq = sum(biomass(sim_start).species[aq_prod])
                    sim_start_bio_int_aq = sum(biomass(sim_start).species[aq_int])
                    sim_start_bio_pred_aq = sum(biomass(sim_start).species[aq_pred])

                # stability
                    sim_start_stability = community_cv(sim_start)

                # shannon
                    sim_start_shannon = shannon_diversity(sim_start)

                # push data to DataFrames
                    push!(biomass_df, [j, "start", "terrestrial", "prod", sim_start_bio_prod_t])
                    push!(biomass_df, [j, "start", "terrestrial", "int", sim_start_bio_int_t])
                    push!(biomass_df, [j, "start", "terrestrial", "pred", sim_start_bio_pred_t])
                    push!(biomass_df, [j, "start", "aquatic", "prod", sim_start_bio_prod_aq])
                    push!(biomass_df, [j, "start", "aquatic", "int", sim_start_bio_int_aq])
                    push!(biomass_df, [j, "start", "aquatic", "pred", sim_start_bio_pred_aq])
                    
                    push!(stability_df, [j, "start", sim_start_stability])
                    
                    push!(shannon_df, [j, "start", sim_start_shannon])

        #### INTRODUCE BEAVER

            # Note: beaver biomass currently discrete value
            # Also, what to do if beaver dies...

            ## Get final biomass values for sim_start
                Bbeaver = sim_start.u[end,:][1]
            
            ## Add beaver biomass to Bbeaver
                Bbeaver[beaver] = 0.25

        #### RUN INTRODUCED COMMUNITY

            ## Run introduced community
                sim_beaver = simulate(params, Bbeaver)

            ## Collect data

                # Biomass
                    sim_beaver_bio_prod_t = sum(biomass(sim_beaver).species[t_prod])
                    sim_beaver_bio_int_t = sum(biomass(sim_beaver).species[t_int])
                    sim_beaver_bio_pred_t = sum(biomass(sim_beaver).species[t_pred])
                    sim_beaver_bio_prod_aq = sum(biomass(sim_beaver).species[aq_prod])
                    sim_beaver_bio_int_aq = sum(biomass(sim_beaver).species[aq_int])
                    sim_beaver_bio_pred_aq = sum(biomass(sim_beaver).species[aq_pred])

                # stability
                    sim_beaver_stability = community_cv(sim_beaver)

                # shannon
                    sim_beaver_shannon = shannon_diversity(sim_beaver)

                # push data to DataFrames
                    push!(biomass_df, [j, "introduced", "terrestrial", "prod", sim_beaver_bio_prod_t])
                    push!(biomass_df, [j, "introduced", "terrestrial", "int", sim_beaver_bio_int_t])
                    push!(biomass_df, [j, "introduced", "terrestrial", "pred", sim_beaver_bio_pred_t])
                    push!(biomass_df, [j, "introduced", "aquatic", "prod", sim_beaver_bio_prod_aq])
                    push!(biomass_df, [j, "introduced", "aquatic", "int", sim_beaver_bio_int_aq])
                    push!(biomass_df, [j, "introduced", "aquatic", "pred", sim_beaver_bio_pred_aq])
                    
                    push!(stability_df, [j, "introduced", sim_beaver_stability])

                    push!(shannon_df, [j, "introduced", sim_beaver_shannon])



        #### START AQUATIC SUCCESSION LOOP
            
            # Note: dispersal chance is currently still set to just 50:50
            # Unsure if code finds coloniser prey or predators (columns vs rows)
            # Currently samples from every aquatic species. Could alter so it samples from non-producers, then if they dont have biomass then a producer can colonise?
            # Omitted code that informed if species was already present or not for now.
            # coloniser starting biomass currently a DISCRETE value

            ## Set timestep interval between succession events
            T_interval = 20
            
            ## Set how many colonisation events occur
            colonisation_number = 40

            ## Set proportion of metacommunity biomass given to coloniser when they appear
            coloniser_biomass_proportion = 0.25

            ## Set max attempts before loop gives update
            max_attempts = 1e6

            ## Create arrays to collect all data required for the loop
                global prevbiomass = []
                global startbiomass = []
                global coloniser = []
                global simulations = []
            
            ## Add sim_beavers final biomass to prevbiomass to give starting biomass values
                push!(prevbiomass, sim_beaver.u[end, :][1])

            ## Start the loop
                for i in 1:S_a  #colonisation_number

                    # inform iteration number
                    println()
                    @info "Iteration number: $i"

                    # reset # of attempts
                    attempts = 1

                    # pick a potential coloniser
                        i_coloniser = sample(setdiff(aquatic_species, coloniser))
                        #i_coloniser = sample(aquatic_species)

                    # Check suitability (1. Is its prey present in the community? 2. 50:50 chance of dispersal)
                    while attempts <= max_attempts #true
                        # identify prey of potential coloniser
                        coloniser_prey = fw_postlinks_matrix[1:length(metacom.species), i_coloniser]
                        coloniser_prey = findall(x -> x == 1, coloniser_prey)

                        # check if the potential coloniser is a producer. If so, then exempt it from needing prey biomass
                        is_prod = i_coloniser in aq_prod

                        # Sum biomass of coloniser's prey
                        coloniser_prey_biomass = sum(prevbiomass[end][coloniser_prey])

                        # check conditions
                        if (is_prod || coloniser_prey_biomass > 0) && rand() > 0.5
                            break
                        else
                            i_coloniser = sample(setdiff(aquatic_species, coloniser))
                            #i_coloniser = sample(aquatic_species)
                            attempts += 1
                        end
                    end

                    # Check if a suitable colonizer was found within the maximum attempts
                    if attempts > max_attempts
                        @warn "No suitable colonizer found after $max_attempts attempts. Exiting the loop."
                        break
                    end

                    # push! successful coloniser into array
                    push!(coloniser, i_coloniser)

                    # Say which species was chosen
                    @info "coloniser species: $(coloniser[i])"

                    functional_group = if i_coloniser in t_prod
                        "terrestrial producer"
                    elseif i_coloniser in t_int
                        "terrestrial intermediate consumer"
                    elseif i_coloniser in t_pred
                        "terrestrial predator"
                    elseif i_coloniser in aq_prod
                        "aquatic producer"
                    elseif i_coloniser in aq_int
                        "aquatic intermediate consumer"
                    else
                        "aquatic predator"
                    end

                    @info "Coloniser species: $(coloniser[i]), $functional_group"

                    # Prepare starting biomasses for the simulation
                    lastbiomass = deepcopy(prevbiomass[end])
                    push!(startbiomass, lastbiomass)

                    # Give the coloniser some mass
                    startbiomass[end][i_coloniser] = 1 #prevbiomass[end][i_coloniser] + (coloniser_biomass_proportion * B_reference[coloniser[end]])

                    # Start simulation
                    sim_i = simulate(params, startbiomass[end], tmax = T_interval, verbose = false)

                    # push the simulation into the simulations array
                    push!(simulations, sim_i)

                    # add the final biomasses to the prevbiomass array for the next iteration
                    push!(prevbiomass, sim_i.u[end, :][1])

                    # Collect data

                        # Biomass
                        sim_i_bio_prod_t = sum(biomass(sim_i).species[t_prod])
                        sim_i_bio_int_t = sum(biomass(sim_i).species[t_int])
                        sim_i_bio_pred_t = sum(biomass(sim_i).species[t_pred])
                        sim_i_bio_prod_aq = sum(biomass(sim_i).species[aq_prod])
                        sim_i_bio_int_aq = sum(biomass(sim_i).species[aq_int])
                        sim_i_bio_pred_aq = sum(biomass(sim_i).species[aq_pred])

                        # stability
                        sim_beaver_stability = community_cv(sim_i)

                        # shannon
                        sim_beaver_shannon = shannon_diversity(sim_i)

                        # push data to DataFrames
                        push!(biomass_df, [j, i, "terrestrial", "prod", sim_beaver_bio_prod_t])
                        push!(biomass_df, [j, i, "terrestrial", "int", sim_beaver_bio_int_t])
                        push!(biomass_df, [j, i, "terrestrial", "pred", sim_beaver_bio_pred_t])
                        push!(biomass_df, [j, i, "aquatic", "prod", sim_beaver_bio_prod_aq])
                        push!(biomass_df, [j, i, "aquatic", "int", sim_beaver_bio_int_aq])
                        push!(biomass_df, [j, i, "aquatic", "pred", sim_beaver_bio_pred_aq])
                        
                        push!(stability_df, [j, i, sim_beaver_stability])

                        push!(shannon_df, [j, i, sim_beaver_shannon])
            
                end
    end

#### EXPORT DATA TO CSV

using CSV

    # Define the folder path
    folder_path = "H:\\.shortcut-targets-by-id\\1elqgJ1X2yIhECbPAwjaVd7hicbF_oj3z\\Will Woof\\Data"

    # Export biomass_df to CSV in the Data folder
    CSV.write(joinpath(folder_path, "biomass_data.csv"), biomass_df)

    # Export stability_df to CSV in the Data folder
    CSV.write(joinpath(folder_path, "stability_data.csv"), stability_df)

    # Export shannon_df to CSV in the Data folder
    CSV.write(joinpath(folder_path, "shannon_data.csv"), shannon_df)






    plot(simulations[30])


    get_extinct_species(simulations[8])
    aq_pred




    coloniser_prey = fw_postlinks_matrix[1:length(metacom.species), metaspecies[19]]
    coloniser_prey = findall(x -> x == 1, coloniser_prey)
    coloniser_prey_biomass = sum(prevbiomass[6][coloniser_prey])
