# This file is part of Kpax3. License is MIT.

function plotk(k::Vector{Int},
               pk::Vector{Float64};
               xlim::UnitRange{Int}=1:0,
               xticks::Vector{Int}=zeros(Int, 0))
  x = if length(xlim) > 0
    collect(xlim)
  else
    collect(max(1, k[1] - 5):(k[end] + 5))
  end

  u = length(k)
  v = length(x)

  y = zeros(Float64, v)

  i = 1
  j = 1
  while (i <= u) && (j <= v)
    if k[i] == x[j]
      y[j] = pk[i]
      i += 1
    end

    j += 1
  end

  ticks = if length(xticks) == 0
    Guide.xticks(ticks=:auto)
  else
    Guide.xticks(ticks=xticks)
  end

  plot(x=x, y=y, Geom.bar,
       Theme(default_color=colorant"black", bar_spacing=1mm, plot_padding=0mm),
       Coord.cartesian(xmin=x[1] - 0.5, xmax=x[end] + 0.5, ymin=0),
       Guide.xlabel("k"),
       Guide.ylabel("p(k | x)"),
       Guide.title(""),
       ticks)
end
