# preference

using EcologicalNetworksDynamics, Plots, Random
Random.seed!(123)

S = 20

# make a food web
# richness, connectance and PPMR
#foodweb = FoodWeb(nichemodel, S; C = 0.1, Z = 10)
foodweb = FoodWeb(4=>(2,1), 3=>(2,1))

# set K and producer competition details
producer_growth_competition = LogisticGrowth(foodweb;
    K = 10,
    a = (diag = 1, offdiag = 1))

# define the functional response and parameters
fr = BioenergeticResponse(foodweb; h = 1.2)

# set up the model
params = ModelParameters(foodweb;
    functional_response = fr,
    producer_growth = producer_growth_competition)

# define initial biomass
B0 = rand(S) # Initial biomass.

# simulate
solution = simulate(params, B0)

# find extinctions
get_extinct_species(solution)

# plot/visualise
plot(solution)
