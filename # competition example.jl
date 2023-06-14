# competition

using EcologicalNetworksDynamics, Plots, Random

Random.seed!(123)

S = 20
foodweb = FoodWeb(nichemodel(S, 0.05))

producer_competition = ProducerCompetition(foodweb; αii = 1, αij = 0.8)
environment = Environment(foodweb, K = 10)
params = ModelParameters(foodweb; producer_competition)
B0 = rand(S) # Initial biomass.
solution = simulate(params, B0; extinction_threshold, tmax, verbose)