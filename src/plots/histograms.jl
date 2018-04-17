# This file is part of Kpax3. License is MIT.

function plotk(k::Vector{Int},
               pk::Vector{Float64};
               xlim::UnitRange{Int}=1:0,
               xticks::Vector{Int}=zeros(Int, 0),
               width::Real=640.0,
               height::Real=360.0)
  x = if length(xlim) > 0
    collect(xlim)
  else
    collect(max(1, k[1] - 5):(k[end] + 5))
  end

  u = length(k)
  v = length(x)

  y = zeros(Float64, v)
  M = 0.0

  i = 1
  j = 1
  while (i <= u) && (j <= v)
    if k[i] == x[j]
      y[j] = pk[i]

      if y[j] > M
        M = y[j]
      end

      i += 1
    end

    j += 1
  end

  multiplier = 10.0
  while (M * multiplier) <= 1.0
    multiplier *= 10
  end
  multiplier *= 10

  M = ceil(M * multiplier) / multiplier

  Plots.bar(x, y,
            # bar options
            fillcolor=:black, orientation=:vertical,
            # plot options
            background_color=:white, html_output_format=:svg,
            size=(width, height), window_title="",
            # subplot options
            grid=true, legend=:none, title="",
            # x axis
            xlabel="k", xlims=(x[1] - 0.5, x[end] + 0.5),
            xticks=length(xticks) == 0 ? :auto : xticks,
            # y axis
            ylabel="p(k | x)", ylims=(0, M))
end
