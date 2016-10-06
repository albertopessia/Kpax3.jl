# This file is part of Kpax3. License is MIT.

function reorderunits(R::Vector{Int},
                      P::Matrix{Float64},
                      clusterorder::Vector{Int})
  n = length(R)
  k = maximum(R)

  M = zeros(Float64, n)
  v = zeros(Float64, n)
  for j in 1:(n - 1), i in (j + 1):n
    if R[i] == R[j]
      M[i] += P[i, j]
      v[i] += 1

      M[j] += P[i, j]
      v[j] += 1
    end
  end

  M ./= v

  neworder = zeros(Int, n)
  midpoint = zeros(Float64, k)
  seppoint = zeros(Float64, k)

  h = 1
  u = 1
  for g in 1:k
    idx = find(R .== clusterorder[g])
    u = length(idx)
    ord = sortperm(M[idx], rev=true)

    copy!(neworder, h, idx[ord], 1, u)
    midpoint[g] = (2 * h + u - 1) / 2
    seppoint[g] = h - 0.5

    h += u
  end

  (neworder, midpoint, seppoint)
end

function plotP(R::Vector{Int},
               P::Matrix{Float64};
               clusterorder::Vector{Int}=zeros(Int, 0),
               clusterlabel::Vector{String}=fill("", 0),
               linesep::Bool=true)
  n = length(R)
  k = maximum(R)

  if length(clusterorder) == 0 || length(clusterorder) != k
    clusterorder = collect(1:k)
  end

  if length(clusterlabel) == 0 || length(clusterlabel) != k
    clusterlabel = map(l -> string(l), clusterorder)
  end

  (ord, mid, sep) = reorderunits(R, P, clusterorder)

  N = n^2

  x = zeros(Int, N)
  y = zeros(Int, N)
  z = fill("", N)

  idx = 1
  pij = 0.0
  for j = 1:n, i = 1:n
    x[idx] = j
    y[idx] = i

    pij = P[ord[i], ord[j]]
    z[idx] = if pij <= 0.2
      "[0.0, 0.2]"
    elseif pij <= 0.4
      "(0.2, 0.4]"
    elseif pij <= 0.6
      "(0.4, 0.6]"
    elseif pij <= 0.8
      "(0.6, 0.8]"
    else
      "(0.8, 1.0]"
    end

    idx += 1
  end

  plot(layer(yintercept=sep, Geom.hline(color=colorant"black")),
       layer(xintercept=sep, Geom.vline(color=colorant"black")),
       layer(x=x, y=y, color=z, Geom.rectbin),
       Coord.cartesian(yflip=true, fixed=true, xmin=0.5, xmax=n + 0.5, ymin=0.5,
                       ymax=n + 0.5),
       Scale.color_discrete_manual(colorant"#FFFFFF", colorant"#CCCCCC",
                                   colorant"#999999", colorant"#666666",
                                   colorant"#333333", colorant"#000000",
                                   levels=["[0.0, 0.2]", "(0.2, 0.4]",
                                           "(0.4, 0.6]", "(0.6, 0.8]",
                                           "(0.8, 1.0]"]),
       Scale.x_continuous(labels=l -> string(clusterlabel[mid .== l][1])),
       Scale.y_continuous(labels=l -> string(clusterlabel[mid .== l][1])),
       Theme(background_color=colorant"#FFFFFF", panel_fill=colorant"#FFFFFF",
             plot_padding=0mm, grid_line_width=0mm,
             line_style=get_stroke_vector(:dash),
             line_width=if linesep 0.3mm else 0mm end),
       Guide.xticks(ticks=mid),
       Guide.yticks(ticks=mid),
       Guide.xlabel(""),
       Guide.ylabel(""),
       Guide.title(""),
       Guide.colorkey(""))
end
