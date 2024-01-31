# Reproduce Fig. 2 of Binzer 2016: https://doi.org/10.1111/gcb.13086
# Ismael update to APB's and Hana's code.
# Using dev version of EcologicalNetworksDynamics.

import Random: seed!
using CairoMakie
using DataFrames
using Distributions
using EcologicalNetworksDynamics
using Statistics

n_cells = 20
T_values = 273.15 .+ LinRange(0, 40, n_cells)
K_values = LinRange(1, 20, n_cells)
n_fw = 5 # Number of food web replicates.
initial_richness = 30
Z = 100 # Predator-prey mass ratio.
C = 0.1

"""
Generate species body masses according to Binzer 2016.
Body masses are generated using the usual scaling of body mass with trophic levels
through the Predator-Prey Mass Ratio (PPMR, `Z`).
In addition, a noise is added on top of the of the trophic levels
"to avoid that all species of a trophic level are equally sized".
"""
function generate_bodymass(foodweb, Z; m0 = 0.01)
    tl = trophic_levels(foodweb)
    noise = rand(Normal(), richness(foodweb))
    m0 * Z .^ (tl .- 1 .+ noise)
end

# Generate food webs.
foodweb_bank = []
for _ in 1:n_fw
    fw = FoodWeb(nichemodel, initial_richness; C)
    fw.M = generate_bodymass(fw, Z)
    push!(foodweb_bank, fw)
end

# Run simulations.
df = DataFrame(;
    eutrophication = Float64[],
    temperature = Float64[],
    persistence = Float64[],
)
tmax = 31536000000 / 100
Threads.@threads for K_idx in 1:length(K_values)
    K0 = K_values[K_idx]
    for (T_idx, T) in enumerate(T_values)
        for (fw_idx, foodweb) in enumerate(foodweb_bank)
            functional_response = ClassicResponse(foodweb; h = 1.2)
            biorates = BioRates(foodweb; d = 0) # Set mortality rate to 0.
            model = ModelParameters(foodweb; functional_response, biorates)
            K = exp_ba_carrying_capacity(; aₚ = K0)
            temperature_scaling = ExponentialBA(; K)
            set_temperature!(model, T, temperature_scaling)
            model_copy = deepcopy(model)
            model_copy.network.M = [1 for _ in 1:initial_richness]
            K_vec = model_copy.producer_growth.K # For consicennes, cf. line below.
            B0_cons = mean(K_vec[producers(foodweb)]) / 8
            B0 = [isproducer(i, foodweb) ? K_vec[i] : B0_cons for i in 1:initial_richness]
            B0 = Float64.(B0)
            callback =
                EcologicalNetworksDynamics.ExtinctionCallback(1e-12, model_copy, false)
            sol = simulate(
                model_copy,
                B0;
                # alg = Tsit5(),
                tmax,
                callback,
                adaptive = true,
                dt = 24 * 60 * 60,
                saveat = 24 * 60 * 60,
            )
            final_richness = richness(sol[end])
            persistence = final_richness / initial_richness
            push!(df, [K0, T, persistence])
        end
        @info "Eutrophication = $K0 - Temperature = $T: done."
    end
end

processed_data = combine(groupby(df, [:eutrophication, :temperature]), :persistence => mean)

# make plot
fig = Figure()
ax = Axis(fig[1, 1])
pp = CairoMakie.heatmap!(ax, processed_data[!,:eutrophication], 
                            processed_data[!,:temperature], 
                            processed_data[!,:persistence_mean])
Colorbar(fig[1,2], pp)
fig

# Plot heatmap.
with_theme() do
    colormap = Reverse(:algae)
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "Eutrophication", ylabel = "Temperature")
    hm = heatmap!(
        ax,
        processed_data.eutrophication,
        processed_data.temperature .- 273.15,
        processed_data.persistence_mean;
        colormap,
        # colorrange = (0, 1),
    )
    # Colorbar(fig[1, 2]; limits = (0, 1), label = "Persistence", colormap)
    #save("/tmp/plot.png", fig)
    hm
end
hm