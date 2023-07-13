using Pkg
using EcologicalNetworksDynamics
using Plots, Random
using Statistics
using LinearAlgebra
using Distributions, DataFrames, StatsPlots, CSV

function async(solution, last = "50%")
        mat = extract_last_timesteps(solution; last, idxs = [3,4])
        cov_mat = cov(mat)
        var_sp = Diagonal(cov_mat)
        var_tot = sum(cov_mat)
        async = sum(sqrt(var_sp)) / sqrt(var_tot)
        async
end

C = range(0.0, 0.15, length = 10) 
B = range(0.0, 0.025, length = 10)  # range from 0.01, 0.001, 0, 0, 0

# Creating a repeat for varying prey preference on 1 consumer
# The kuramoto food web 
foodweb = FoodWeb([4=>(2, 1), 3=>(2, 1)])

df_collect_comb = DataFrame(B = [], C = [], CoefficientofVariation = [], LoreauSync = [], ReumanSync = [])


for i in 1:length(C)
        for j in 1:length(B)

                # BIORATES
                # set mortality of consumers and producer growth rate 
                # Kuramoto values
                br = BioRates(foodweb,
                # consumers have death rate
                d = [0, 0, 0.1, 0.1],
                # producers have growth rates 
                r = [0.3, 0.3, 0, 0],
                # efficiencies 100%
                e = ones(4,4),
                # metabolism and max cons = 0 (?)
                x = 0.0,
                y=0.0)

                # FUNCTIONAL RESPONSE and PREFERENCE
                fr = ClassicResponse(foodweb;
                # hill exponent
                h = 1,
                # handling and attack
                hₜ = 3, aᵣ = 0.7,
                # consumer interference
                c = 0.0,
                # consumer preference
                ω = [0 0 0 0;
                0 0 0 0;
                1.0 B[j] 0 0;
                B[j] 1.0 0 0])

                # PRODUCER GROWTH, CARRYING CAPACITY AND COMPETITION
                producer_growth_competition = LogisticGrowth(foodweb; 
                        a = (diag = 1, offdiag = C[i]), K = 1)

                # THE MODEL 
                params = ModelParameters(foodweb;
                        biorates = br,
                        functional_response = fr,
                        producer_growth = producer_growth_competition)

                ###### MODELLING ######

                # define initial biomass (Kuramoto starters)
                B0 = [0.1, 0.3, 0.1, 0.3] # Initial biomass

                # simulate
                solution = simulate(params, B0, tmax = 600)
                
                # collect data 
                
                Rsync = async(solution, "25%")
                var= coefficient_of_variation(solution).community
                Lsync = coefficient_of_variation(solution, idxs = [3,4], last = "25%").synchrony
                # write to df        
                push!(df_collect_comb, [B[j], C[i], var, Lsync, Rsync])
        end
end

df_collect_comb

# # to use in R
# CSV.write("karimoto_basics.csv", df_collect_comb)

# subset to B = 0 for C effect and C = 0 for B effect
C_only = filter(:B => ==(0), df_collect_comb)
B_only = filter(:C => ==(0), df_collect_comb)

p1 = @df C_only plot(:C, :LoreauSync, xlabel = "C", title = "Loreau")
p2 = @df C_only plot(:C, :ReumanSync, xlabel = "C", title  = "Reuman")
p3 = @df B_only plot(:B, :LoreauSync, xlabel = "B")
p4 = @df B_only plot(:B, :ReumanSync, xlabel = "B")

# combined plot
plot(p1, p2, p3, p4, ncol = 2)