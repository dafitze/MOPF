using DataFrames
using Turing

using GLMakie, AlgebraOfGraphics
using AlgebraOfGraphics: density

# ===============================================
# Simulated Data
# ===============================================
include("data_simulation.jl")
include("plot_pf.jl")

d_sim = simulate_data()
plot_pf(d_sim)

# model
# -----------------------------------------------
@model function logreg_one(y,x)
  # prior
  α ~ Normal(0,10)
  β₁ ~ Normal(0,10)
  linpred = α .+ β₁ * x
  # likelihood
  y .~ Bernoulli.(logistic.(linpred))
end


# prior
# -----------------------------------------------
prior_fit = sample(logreg_one(d_sim.response, d_sim.stimulus),
                   Prior(),
                   2000)

#model_missing = logreg_one(similar(d_sim.response,Missing), d_sim.stimulus) # instantiate the "predictive model
#prior_check = predict(model_missing, prior_fit);

# posterior fit
# -----------------------------------------------
posterior_fit = sample(logreg_one(d_sim.response, d_sim.stimulus),
                       NUTS(),
                       MCMCThreads(),
                       2000,
                       4)



draw(p * mapping())

  pars = names(posterior_fit, :parameters)
chain_mapping =
  mapping(pars .=> "sample value") *
  mapping(; color=:chain => nonnumeric, row=dims(1) => renamer(pars))


d1 = posterior_fit
d2 = pars_sim
layers1 = data(posterior_fit) * density()
layers2 = data(df2) * visual(Scatter)

draw(layers1 * chain_mapping)

draw(layers2 * chain_mapping)

df1 = (x=rand(100), y=rand(100), i=rand(["a", "b", "c"], 100))
df2 = (x=[0, 1], y=[0.5, 0.5], i=fill("b", 2))
layers1 = data(df1) * visual(Scatter) 
layers2 = data(df2) * visual(Lines)
draw(layers1 * mapping(:x, :y, row = :i) + layers2 * mapping(:x, :y, row=:i))


function plot_pars(fit)
  params = names(fit, :parameters)
  chain_mapping =
  mapping(params .=> "sample value") *
  mapping(; color=:chain => nonnumeric, row=dims(1) => renamer(params))
  p_fit = data(fit) * chain_mapping * density()
  f = Figure(; resolution=(800, 600))
  draw!(f[1, 1], p_fit; axis=(; ylabel="density"))
  return p_fit
end

p = plot_pars(posterior_fit)
draw(p)

pars_sim = DataFrame(pars = ["β₀", "β₁"],
                     xs =  [unique(d_sim.β₀)[1], unique(d_sim.β₁)[1]],
                     ys = [0,0])



f = Figure(; resolution=(800, 600))
ax_param = Axis(f[1,1])
#xlabel = "",
#ylabel = "Density")
draw!(f[1,1], p)

function plot_param(fit, simulated_data)
  #pars_sim = DataFrame(pars = ["β₀", "β₁"],
  #                     val =  [unique(simulated_data.β₀)[1], unique(simulated_data.β₁)[1]])
  pars_fit = names(fit, :parameters)
  #
  chain_mapping =
  mapping(pars_fit .=> "sample value") *
  mapping(; color=:chain => nonnumeric, row=dims(1) => renamer(pars_fit))
  #
  dens = data(fit) * chain_mapping * density()
  #sim = data(pars_sim) * mapping(:val, 1.0) * visual(Scatter, markersize = 50)
  #
  f = Figure(; resolution=(800, 600))
  ax_param = Axis(f[1,1],
                  xlabel = "",
                  ylabel = "Density")
  draw!(ax_param, dens; facet = (;linkxaxes = :none))
  #draw!(ax_param, sim)#; facet = (;linkxaxes = :none))
  return f
end
