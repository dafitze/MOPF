using Chain, DataFrames
using StatsFuns: logistic
using Distributions

function simulate_data(
    β₀ = 0.0,
    β₁ = 2.0,
    γ = 0,
    λ = 0,
    vpn = 1,
    nreps = 20)


  stimulus = [-214,-180,-146,-112,-78,-44,-10,10,44,78,112,146,180,214]./100

  df = @chain Iterators.product(vpn,1:nreps, stimulus, β₀, β₁, γ, λ) begin
    collect()
    reduce(vcat, _)
    DataFrame([:vpn, :rep, :stimulus, :β₀, :β₁, :γ, :λ])
    DataFrames.transform(_, [:β₀,:β₁,:γ,:λ,:stimulus] => ByRow((b0,b1,guess,lapse,stim) -> (1 - guess - lapse) * logistic(b0+(b1*stim))) => :θ)
    DataFrames.transform(_, :θ => ByRow(p -> rand(Binomial(1,p))) => :response)
    sort(_, [:vpn,:stimulus,:rep])
  end
  return df
end
