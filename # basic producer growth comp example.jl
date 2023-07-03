# competition

using EcologicalNetworksDynamics, Plots, Random

Random.seed!(123)

S = 20

foodweb = FoodWeb(nichemodel(S, 0.05))

producer_growth_competition = LogisticGrowth(foodweb; K = 10, a = (diag = 1, offdiag = 0.8))
fr = BioenergeticResponse(foodweb; h = 1.2)

params = ModelParameters(foodweb;
    functional_response = fr,
    producer_growth = producer_growth_competition)

B0 = rand(S) # Initial biomass.
solution = simulate(params, B0)

get_extinct_species(solution)

plot(solution)
