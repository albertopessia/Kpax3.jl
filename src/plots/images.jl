# This file is part of Kpax3. License is MIT.

function plotP(R::Vector{Int},
               P::Matrix{Float64};
               clusterorder::Vector{Int}=zeros(Int, 0),
               clusterlabel::Vector{String}=fill("", 0),
               linesep::Bool=true,
               width::Real=640.0,
               height::Real=640.0)
  n = length(R)
  k = maximum(R)

  if length(clusterorder) == 0 || length(clusterorder) != k
    clusterorder = collect(1:k)
  end

  if length(clusterlabel) == 0 || length(clusterlabel) != k
    clusterlabel = map(l -> string(l), clusterorder)
  end

  (ord, mid, sep) = reorderunits(R, P, clusterorder)

  mid .-= 0.5
  sep .-= 0.5

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
      "#FFFFFF"
    elseif pij <= 0.4
      "#F0F0F0"
    elseif pij <= 0.6
      "#BDBDBD"
    elseif pij <= 0.8
      "#636363"
    else
      "#000000"
    end

    idx += 1
  end

  # try to minimize the number of rectangles. We will do this greedly by
  # dividing the image into maximal adjacent rectangles
  # Hopefully the SVG file won't be too big
  #
  # TODO: improve the algorithm
  processed = falses(n, n)

  xmin = zeros(Int, 0)
  xmax = zeros(Int, 0)
  ymin = zeros(Int, 0)
  ymax = zeros(Int, 0)
  col = fill("", 0)

  # elements on the diagonal have always p = 1.0
  val = "#000000"

  j = 1
  while j <= n
    jmax = expandsquarediag(j, n, val, z)

    push!(xmin, j - 1)
    push!(xmax, jmax)
    push!(ymin, j - 1)
    push!(ymax, jmax)
    push!(col, val)

    processed[j:jmax, j:jmax] .= true
    j = jmax + 1
  end

  # scan the lower triangular matrix
  j = 1
  while j < n
    # start a new column from the first non processed element
    i = j + 1
    while i <= n
      if !processed[i, j]
        val = z[LinearIndices((n, n))[i, j]]

        (imax, jmax) = expandrect(i, j, n, n, val, z, processed)

        # lower triangular
        push!(xmin, j - 1)
        push!(xmax, jmax)
        push!(ymin, i - 1)
        push!(ymax, imax)
        push!(col, val)

        # upper triangular (transpose)
        push!(xmin, i - 1)
        push!(xmax, imax)
        push!(ymin, j - 1)
        push!(ymax, jmax)
        push!(col, val)

        processed[i:imax, j:jmax] .= true
        i = imax + 1
      else
        while i <= n && processed[i, j]
          i += 1
        end
      end
    end

    j += 1
  end

  xborder = zeros(Float64, 5, k)
  yborder = zeros(Float64, 5, k)
  for i = 1:k
    xborder[:, i] = [sep[i]; sep[i]; sep[i + 1]; sep[i + 1]; sep[i]]
    yborder[:, i] = [sep[i]; sep[i + 1]; sep[i + 1]; sep[i]; sep[i]]
  end

  rectangles = reshape([Plots.Shape([(xmin[i], ymin[i]),
                                     (xmin[i], ymax[i]),
                                     (xmax[i], ymax[i]),
                                     (xmax[i], ymin[i])])
                        for i = 1:length(xmin)], 1, length(xmin))

  pcol = reshape(col, 1, length(col))
  # create an empty plot and then add all the rectangles on top of it
  mcol = ["#FFFFFF" "#F0F0F0" "#BDBDBD" "#636363" "#000000"]
  mlab = ["[0.0, 0.2]" "(0.2, 0.4]" "(0.4, 0.6]" "(0.6, 0.8]" "(0.8, 1.0]"]

  p = Plots.plot(xborder, yborder,
                 # cluster borders
                 seriestype=:path, linecolor=:black,
                 # plot options
                 background_color=:white, html_output_format=:svg,
                 size=(width, height), window_title="",
                 # subplot options
                 grid=false, label="", legend=:none, title="",
                 # axes
                 formatter=l -> string(clusterlabel[mid .== l][1]),
                 # x axis
                 xlabel="Samples by cluster", xlims=(0, n), xticks=mid,
                 # y axis
                 ylabel="Samples by cluster", ylims=(0, n), yticks=mid,
                 yflip=true)
  Plots.plot!(p, fill(-2, 2, length(mcol)), fill(-0.5, 2, length(mcol)),
              seriestype=:scatter, markershape=:rect, markercolor=mcol,
              label=mlab, legend=:right)
  Plots.plot!(p, rectangles, fillcolor=pcol, linealpha=0.0, linewidth=0.0,
              linestyle=:dot, label="")

  p
end

function plotC(site::Vector{Int},
               freq::Vector{Float64},
               C::Matrix{Float64};
               width::Real=640.0,
               height::Real=360.0)
  m = site[end]
  M = length(site)

  noise = zeros(Float64, m)
  weak = zeros(Float64, m)

  key = 1
  for b in 1:M
    if site[b] != key
      weak[key] += noise[key]
      key += 1
    end

    noise[key] += freq[b] * C[1, b]
    weak[key] += freq[b] * C[2, b]
  end

  weak[key] += noise[key]

  mcol = ["#FFFFFF" "#808080" "#000000"]
  mlab = ["Noise" "Weak signal" "Strong signal"]

  p = Plots.plot(1:m, ones(Float64, m),
                 seriestype=:path, linecolor=:black, fill=(0, :black),
                 # plot options
                 background_color=:white, html_output_format=:svg,
                 size=(width, height), window_title="",
                 # subplot options
                 grid=false, label="", legend=:none, title="",
                 # x axis
                 xlabel="Site", xlims=(1, m),
                 # y axis
                 ylabel="Posterior probability", ylims=(0, 1))
  Plots.plot!(p, 1:m, weak,
              seriestype=:path,
              linecolor="#808080",
              fill=(0, "#808080"),
              label="", legend=:none)
  Plots.plot!(p, 1:m, noise,
              seriestype=:path, linecolor=:white, fill=(0, :white),
              label="", legend=:none)
  Plots.plot!(p, fill(-2, 2, length(mcol)), fill(-0.5, 2, 3),
              seriestype=:scatter, markershape=:rect, markercolor=mcol,
              label=mlab, legend=:right)

  p
end

function plotD(x::AminoAcidData,
               state::AminoAcidState;
               clusterorder::Vector{Int}=zeros(Int, 0),
               clusterlabel::Vector{String}=fill("", 0),
               linesep::Bool=true,
               width::Real=640.0,
               height::Real=360.0)
  (m, n) = size(x.data)
  M = length(x.ref)

  if length(clusterorder) == 0 || length(clusterorder) != state.k
    clusterorder = collect(1:state.k)
  end

  if length(clusterlabel) == 0 || length(clusterlabel) != state.k
    clusterlabel = map(l -> string(l), clusterorder)
  end

  v = zeros(Float64, state.k)
  for i in 1:n
    v[state.R[i]] += 1.0
  end
  v /= n

  v = [0.0; v[clusterorder]]

  for g in 2:(state.k + 1)
    v[g] += v[g - 1]
  end

  xax = collect(1:M)

  yax = zeros(Float64, state.k)
  for g in 1:state.k
    yax[g] = (v[g] + v[g + 1]) / 2
  end

  # colors
  colset = [("A", "Ala", "#FF3232");
            ("R", "Arg", "#FF7F00");
            ("N", "Asn", "#CAB2D6");
            ("D", "Asp", "#FDBF6F");
            ("C", "Cys", "#33A02C");
            ("E", "Glu", "#1F4E78");
            ("Q", "Gln", "#FF3300");
            ("G", "Gly", "#660066");
            ("H", "His", "#B2DF8A");
            ("I", "Ile", "#CC6600");
            ("L", "Leu", "#375623");
            ("K", "Lys", "#A6CEE3");
            ("M", "Met", "#CC3399");
            ("F", "Phe", "#FB9A99");
            ("P", "Pro", "#C81414");
            ("S", "Ser", "#1F78B4");
            ("T", "Thr", "#6A3D9A");
            ("W", "Trp", "#FFFF99");
            ("Y", "Tyr", "#FF66CC");
            ("V", "Val", "#993300");
            ("U", "Sec", "#00FFFF");
            ("O", "Pyl", "#FFFF00");
            ("B", "Asx", "#E4B9A3");
            ("Z", "Glx", "#8F413C");
            ("J", "Xle", "#825E12");
            ("X", "Xaa", "#000000");
            ("-",   "-", "#E0E0E0");
            ("*",   "*", "#000000");
            ( "",    "", "#FFFFFF")]

  coltable = Dict(map(a -> (a[1], a[3]), colset))

  # indices
  idx = 1
  s = 0
  t = 1
  b = 1
  g = 1
  h = 1

  # temporary values
  c = 0x01
  w = 0
  flag = false

  z = fill("#FFFFFF", state.k * M)

  for j in 1:M
    if x.ref[j] == UInt8('.')
      s += 1
      c = 0x01

      while (b <= m) && (x.key[b] == s)
        if state.C[state.cl[1], b] > c
          c = state.C[state.cl[1], b]
        end
        b += 1
      end

      if c > 0x02
        for g in 1:state.k
          h = clusterorder[g]
          w = t
          flag = false

          while !flag && w < b
            flag = state.C[state.cl[h], w] == 0x04
            w += 1
          end

          if flag
            z[idx] = coltable[uppercase(string(Char(x.val[w - 1])))]
          end

          idx += 1
        end
      else
        idx += state.k
      end

      t = b
    else
      idx += state.k
    end
  end

  # try to minimize the number of rectangles. We will do this greedly by
  # dividing the image into maximal adjacent rectangles
  # Hopefully the SVG file won't be too big
  #
  # TODO: improve the algorithm
  processed = falses(state.k, M)

  xmin = zeros(Int, 0)
  xmax = zeros(Int, 0)
  ymin = zeros(Float64, 0)
  ymax = zeros(Float64, 0)
  col = fill("", 0)

  j = 1
  while j <= M
    # start a new column from the first non processed element
    i = 1
    while i <= state.k
      if !processed[i, j]
        val = z[LinearIndices((state.k, M))[i, j]]

        (imax, jmax) = expandrect(i, j, state.k, M, val, z, processed)

        push!(xmin, j - 1)
        push!(xmax, jmax)
        push!(ymin, v[i])
        push!(ymax, v[imax + 1])
        push!(col, val)

        processed[i:imax, j:jmax] .= true
        i = imax + 1
      else
        while i <= state.k && processed[i, j]
          i += 1
        end
      end
    end

    j += 1
  end

  rectangles = reshape([Plots.Shape([(xmin[i], ymin[i]),
                                     (xmin[i], ymax[i]),
                                     (xmax[i], ymax[i]),
                                     (xmax[i], ymin[i])])
                        for i = 1:length(xmin)], 1, length(xmin))

  pcol = reshape(col, 1, length(col))
  mcol = reshape(map(a -> a[3], colset[1:26]), 1, 26)
  mlab = reshape(map(a -> a[2], colset[1:26]), 1, 26)

  p = Plots.plot(rectangles,
                 fillcolor=pcol, linecolor=pcol,
                 # plot options
                 background_color=:white, html_output_format=:svg,
                 size=(width, height), window_title="",
                 # subplot options
                 grid=false, label="", legend=:none, title="",
                 # x axis
                 xlabel="Site", xlims=(1, M),
                 # y axis
                 ylabel="Cluster", ylims=(0, 1), yticks=yax, yflip=true,
                 yformatter=l -> string(clusterlabel[yax .== l][1]))
  Plots.plot!(p, fill(-2, 2, length(mcol)), fill(-0.5, 2, length(mcol)),
              seriestype=:scatter, markershape=:rect, markercolor=mcol,
              label=mlab, legend=:right)
  Plots.plot!(p, v[2:state.k], seriestype=:hline, linecolor=:black, label="")

  p
end
