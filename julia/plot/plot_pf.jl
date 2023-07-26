function plot_pf(dat)
  set_aog_theme!()
  update_theme!(fontsize = 40)
  f = Figure(resolution = (1600, 1600))
  ax_pf = Axis(f[1,1],
               # title = "lkj",
               xlabel = "Stimulus",
               ylabel = "P(right|stimulus)",
               yticks = 0:0.5:1)

  ϕ(x, β₀, β₁) = logistic(β₀ + (β₁ * x))

  xs = range(minimum(dat.stimulus), maximum(dat.stimulus), length=251)
  ys = ϕ.(xs, unique(dat.β₀), unique(dat.β₁))
  pf = data((x=xs, y=ys)) * visual(Lines, linewidth = 5) * mapping(:x, :y)
  d_summary = @chain dat begin
    groupby(:stimulus)
    combine(:response => mean => :mean_response)
  end
  p_summary = data(d_summary) * visual(Scatter, markersize = 30, color = :red) * mapping(:stimulus, :mean_response)

  draw!(ax_pf, pf)
  draw!(ax_pf, p_summary)

  draw(pf+p_summary)
  return f
end


