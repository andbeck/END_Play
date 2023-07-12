# competition

using EcologicalNetworksDynamics, Plots, Random

Random.seed!(123)

# The kuramoto food web
foodweb = FoodWeb([4=>(2, 1), 3=>(2, 1)])

# set mortality of consumers and and producer growth rate
# Kuramoto values
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


######## Set Couplings #########

# B = omega value/prefernce/interaction strengths
# C = interspecific competition

# Settings for NO coupling
C = 0 # range from 0 (no plant comp) - 1.2 (inter>intra)
B = 0 # no consumer coupling - the consumers are specialist

# Settings for Resource Coupling
C = 0.2 # range from 0 (no plant comp) - 1.2 (inter>intra)
B = 0 # no consumer coupling - the consumers are specialist

# Settings for Consumer Coupling
# note small values of B => highly asymetric interaction strength.
C = 0
B = 0.01 # range from 0 (no resource coupling) - 0.1 (90:10 ratio)


# Set details of Classic FR following Kuramoto
# hill = 1; ht and attack specific, consumer interference = 0
# note preference matrix with B value
# 0 means no coupling by consumer (they are specialists)
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
# diag is intraspecific and offdiag is interspecific - defines coupling by resources
# massive bias to intraspecific == stabilising.
# start with offdiag = 0.1
producer_growth_competition = LogisticGrowth(foodweb;
    K = 1,
    a = (diag = 1, offdiag = C))


# set up the model
params = ModelParameters(foodweb;
    biorates = br,
    functional_response = fr,
    producer_growth = producer_growth_competition)

    # define initial biomass (Kuramoto starters)
B0 = [0.1, 0.3, 0.1, 0.3] # Initial biomass.

# simulate
solution = simulate(params, B0, tmax = 600)

# plot/visualise consumer behaviour
plot(solution, idxs = [3,4])

# explore synchrony for consumers
coefficient_of_variation(solution, idxs = [3,4]).synchrony

# correct
using Statistics, LinearAlgebra

# Loreau Method
test = transpose(extract_last_timesteps(solution, last = 100, idxs = [3,4]))
cc = sum(cov(test))
ss = sum(std.(eachcol(test)))
cc/(ss^2)

# Reuman Method for Synchrony
mat = transpose(extract_last_timesteps(solution, last = 80, idxs = [3,4]))
cov_mat = cov(mat)
var_sp = Diagonal(cov_mat)
var_tot = sum(cov_mat)
sum(sqrt(var_sp)) / sqrt(var_tot)


#asynchrony = function(mat) {
cov_mat <- cov(mat)
var_sp <- diag(cov_mat)
var_tot <- sum(cov_mat)
​sum(sqrt(var_sp)) / sqrt(var_tot)
#}

function async(solution, last = "50%")
    mat = extract_last_timesteps(solution; last, idxs = [3,4])
    cov_mat = cov(mat)
    var_sp = Diagonal(cov_mat)
    var_tot = sum(cov_mat)
    async = sum(sqrt(var_sp)) / sqrt(var_tot)
    async
end

async(solution)
async(solution, 80)
async(solution, "10%")

# explore community stability
coefficient_of_variation(solution).community

plot(solution)
