using EcologicalNetworksDynamics
using DiffEqCallbacks
using DifferentialEquations

# Disconnected consumer that goes extinct, all fine
A = [0 0 0 0; 0 0 0 0; 1 0 0 0; 0 1 0 0]
foodweb = FoodWeb(A);
params = ModelParameters(foodweb);
B0 = [0.5, 0, .5, .5];

simulate(params, B0)

# Set callback
steady_state_call_back = CallbackSet(TerminateSteadyState(1e-6, 1e-4), ExtinctionCallback(1e-6, params, true),)

tmax_call_back = CallbackSet(
                             ExtinctionCallback(1e-6, params, true),
                            )

# This simulation stops at steady state, i.e. before tmax
sol = simulate(params, B0;
               tmax = 1000,
               callback = steady_state_call_back,
              );
maximum(sol.t)

# This simulation stops at tmax
sol = simulate(params, B0;
               tmax = 1000,
               callback = tmax_call_back,
              )

maximum(sol.t)
