# This file is part of Kpax3. License is MIT.

function plotk(k::Vector{Int},
               pk::Vector{Float64};
               xlim::UnitRange{Int}=1:0,
               xticks::Vector{Int}=zeros(Int, 0),
               width::Real=89.0,
               height::Real=55.0)
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

  ticks = if length(xticks) == 0
    Guide.xticks(ticks=:auto)
  else
    Guide.xticks(ticks=xticks)
  end

  majorfontsize = computefontsize(width, height)
  minorfontsize = 0.8 * majorfontsize
  keytitlefontsize = majorfontsize
  keylabelfontsize = 0.8 * keytitlefontsize

  plot(x=x, y=y, Geom.bar,
       Coord.cartesian(xmin=x[1] - 0.5, xmax=x[end] + 0.5, ymin=0, ymax=M),
       Theme(default_color=colorant"black", background_color=colorant"white",
             panel_fill=colorant"white", bar_spacing=1mm, plot_padding=0mm,
             lowlight_color=c->c, lowlight_opacity=0.0,
             major_label_font_size=majorfontsize,
             minor_label_font_size=minorfontsize,
             key_title_font_size=keytitlefontsize,
             key_label_font_size=keylabelfontsize),
       Guide.xlabel("k", orientation=:horizontal),
       Guide.ylabel("p(k | x)", orientation=:vertical),
       Guide.title(""),
       ticks)
end
