# competition

using EcologicalNetworksDynamics, Plots, Random

Random.seed!(123)

#S = 20

# make a food web
# richness, connectance and PPMR
# foodweb = FoodWeb(nichemodel, S; C = 0.1, Z = 10)

# kuramoto food web
foodweb = FoodWeb([4=>(2, 1), 3=>(2, 1)])

# set mortality of consumers and and producer growth rate
br = BioRates(foodweb,
    # consumers have death rate
    d = [0, 0, 0.1, 0.1],
    # producers have growth rates
    r = [0.3, 0.3, 0, 0],
    # effeciencies 100%
    e = ones(4,4),
    # metabolism and max cons = 0 (?)
    x = 0.0,
    y = 0.0)

# define the functional response and parameters
# set preferences via omega off diagonal

# B = omega value
# Set to 0 for no consumer coupling (e.g. preference)
# 0.01 for very asymetric but coupled
# bigger

B = 0.01

fr = ClassicResponse(foodweb;
    # hill exponent
    h = 1,
    # handling and attack
    hₜ = 3, aᵣ = 0.7,
    # comsumer interference
    c = 0.0,
    # consumer preference
    ω = [0 0 0 0;
    0 0 0 0;
    1.0 B 0 0;
    B 1.0 0 0])

# set K and producer competition details
# diag is intraspecific and offdiag is interspecific
# massive bias to intraspecific == stabilising.
producer_growth_competition = LogisticGrowth(foodweb;
    K = 1,
    a = (diag = 1, offdiag = 0))


# set up the model
params = ModelParameters(foodweb;
    biorates = br,
    functional_response = fr,
    producer_growth = producer_growth_competition)

# define initial biomass
B0 = [0.1, 0.3, 0.1, 0.3] # Initial biomass.

# simulate
solution = simulate(params, B0, tmax = 600)

# find extinctions
get_extinct_species(solution)

# plot/visualise consumer behaviour
plot(solution, idxs = [3,4])

coefficient_of_variation(solution, idxs = [3,4])
