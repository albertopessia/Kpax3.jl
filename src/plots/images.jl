# This file is part of Kpax3. License is MIT.

function plotP(R::Vector{Int},
               P::Matrix{Float64};
               clusterorder::Vector{Int}=zeros(Int, 0),
               clusterlabel::Vector{String}=fill("", 0),
               linesep::Bool=true,
               width::Real=89.0,
               height::Real=89.0)
  n = length(R)
  k = maximum(R)

  if length(clusterorder) == 0 || length(clusterorder) != k
    clusterorder = collect(1:k)
  end

  if length(clusterlabel) == 0 || length(clusterlabel) != k
    clusterlabel = map(l -> string(l), clusterorder)
  end

  (ord, mid, sep) = reorderunits(R, P, clusterorder)

  mid -= 0.5
  sep -= 0.5

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
  val = "(0.8, 1.0]"

  j = 1
  while j <= n
    jmax = expandsquarediag(j, n, val, z)

    push!(xmin, j - 1)
    push!(xmax, jmax)
    push!(ymin, j - 1)
    push!(ymax, jmax)
    push!(col, val)

    processed[j:jmax, j:jmax] = true
    j = jmax + 1
  end

  # scan the lower triangular matrix
  j = 1
  while j < n
    # start a new column from the first non processed element
    i = j + 1
    while i <= n
      if !processed[i, j]
        val = z[sub2ind((n, n), i, j)]

        (imax, jmax) = expand(i, j, n, n, val, z, processed)

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

        processed[i:imax, j:jmax] = true
        i = imax + 1
      else
        while i <= n && processed[i, j]
          i += 1
        end
      end
    end

    j += 1
  end

  majorfontsize = computefontsize(width, height)
  minorfontsize = 0.8 * majorfontsize

  keytitlefontsize = majorfontsize
  keylabelfontsize = 0.8 * keytitlefontsize
  padding = 1.0 * majorfontsize

  Gadfly.plot(Gadfly.layer(x=repeat(sep[1:k], outer=2),
                           y=[sep[1:k]; sep[2:end]],
                           xend=repeat(sep[2:end], outer=2),
                           yend=[sep[1:k]; sep[2:end]],
                           Gadfly.Geom.segment),
              Gadfly.layer(x=[sep[1:k]; sep[2:end]],
                           y=repeat(sep[1:k], outer=2),
                           xend=[sep[1:k]; sep[2:end]],
                           yend=repeat(sep[2:end], outer=2),
                           Gadfly.Geom.segment),
              Gadfly.layer(x_min=xmin, x_max=xmax, y_min=ymin, y_max=ymax,
                           color=col, Gadfly.Geom.rect),
              Gadfly.Coord.cartesian(yflip=true, fixed=true,
                                     xmin=0, xmax=n,
                                     ymin=0, ymax=n),
              Gadfly.Scale.color_discrete_manual(
                Gadfly.@colorant_str("#FFFFFF"),
                Gadfly.@colorant_str("#CCCCCC"),
                Gadfly.@colorant_str("#999999"),
                Gadfly.@colorant_str("#666666"),
                Gadfly.@colorant_str("#333333"),
                Gadfly.@colorant_str("#000000"),
                levels=["[0.0, 0.2]", "(0.2, 0.4]",
                        "(0.4, 0.6]", "(0.6, 0.8]",
                        "(0.8, 1.0]"]),
              Gadfly.Scale.x_continuous(labels=l ->
                string(clusterlabel[mid .== l][1])),
              Gadfly.Scale.y_continuous(labels=l ->
                string(clusterlabel[mid .== l][1])),
              Gadfly.Theme(default_color=Gadfly.@colorant_str("black"),
                           background_color=Gadfly.@colorant_str("white"),
                           panel_fill=Gadfly.@colorant_str("white"),
                           grid_line_width=0Measures.mm,
                           line_width=linesep ? 0.3Measures.mm : 0Measures.mm,
                           lowlight_color=c->c,
                           lowlight_opacity=0.0,
                           major_label_font_size=majorfontsize,
                           minor_label_font_size=minorfontsize,
                           key_title_font_size=keytitlefontsize,
                           key_label_font_size=keylabelfontsize,
                           plot_padding=padding),
              Gadfly.Guide.xticks(ticks=mid),
              Gadfly.Guide.yticks(ticks=mid),
              Gadfly.Guide.xlabel("Samples by cluster"),
              Gadfly.Guide.ylabel("Samples by cluster"),
              Gadfly.Guide.title(""),
              Gadfly.Guide.colorkey(""))
end

function plotC(site::Vector{Int},
               freq::Vector{Float64},
               C::Matrix{Float64};
               width::Real=183.0,
               height::Real=92.0)
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

  majorfontsize = computefontsize(width, height)
  minorfontsize = 0.8 * majorfontsize

  keytitlefontsize = majorfontsize
  keylabelfontsize = 0.8 * keytitlefontsize
  padding = 1.0 * majorfontsize

  Gadfly.plot(Gadfly.layer(x=repeat(1:m, outer=3),
                           ymin=[zeros(Float64, m); noise; weak],
                           ymax=[noise; weak; ones(Float64, m)],
                           color=[fill("Noise", m);
                                  fill("Weak", m);
                                  fill("Strong", m)],
                           Gadfly.Geom.ribbon),
              Gadfly.Coord.cartesian(xmin=1, xmax=m, ymin=0, ymax=1),
              Gadfly.Scale.x_continuous(minvalue=1, maxvalue=m),
              Gadfly.Scale.y_continuous(minvalue=0, maxvalue=1),
              Gadfly.Scale.color_discrete_manual(
                Gadfly.@colorant_str("white"),
                Gadfly.@colorant_str("gray"),
                Gadfly.@colorant_str("black")),
              Gadfly.Theme(default_color=Gadfly.@colorant_str("black"),
                           background_color=Gadfly.@colorant_str("white"),
                           panel_fill=Gadfly.@colorant_str("white"),
                           grid_line_width=0Measures.mm,
                           lowlight_color=c->c,
                           lowlight_opacity=0.0,
                           major_label_font_size=majorfontsize,
                           minor_label_font_size=minorfontsize,
                           key_title_font_size=keytitlefontsize,
                           key_label_font_size=keylabelfontsize,
                           plot_padding=padding),
              Gadfly.Guide.xlabel("Site"),
              Gadfly.Guide.ylabel("Posterior probability"),
              Gadfly.Guide.title(""),
              Gadfly.Guide.colorkey("Signal"))
end

function plotD(x::AminoAcidData,
               state::AminoAcidState;
               clusterorder::Vector{Int}=zeros(Int, 0),
               clusterlabel::Vector{String}=fill("", 0),
               linesep::Bool=true,
               width::Real=183.0,
               height::Real=92.0)
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
  colset = [("A", "Ala", Gadfly.@colorant_str("#FF3232"));
            ("R", "Arg", Gadfly.@colorant_str("#FF7F00"));
            ("N", "Asn", Gadfly.@colorant_str("#CAB2D6"));
            ("D", "Asp", Gadfly.@colorant_str("#FDBF6F"));
            ("C", "Cys", Gadfly.@colorant_str("#33A02C"));
            ("E", "Glu", Gadfly.@colorant_str("#1F4E78"));
            ("Q", "Gln", Gadfly.@colorant_str("#FF3300"));
            ("G", "Gly", Gadfly.@colorant_str("#660066"));
            ("H", "His", Gadfly.@colorant_str("#B2DF8A"));
            ("I", "Ile", Gadfly.@colorant_str("#CC6600"));
            ("L", "Leu", Gadfly.@colorant_str("#375623"));
            ("K", "Lys", Gadfly.@colorant_str("#A6CEE3"));
            ("M", "Met", Gadfly.@colorant_str("#CC3399"));
            ("F", "Phe", Gadfly.@colorant_str("#FB9A99"));
            ("P", "Pro", Gadfly.@colorant_str("#C81414"));
            ("S", "Ser", Gadfly.@colorant_str("#1F78B4"));
            ("T", "Thr", Gadfly.@colorant_str("#6A3D9A"));
            ("W", "Trp", Gadfly.@colorant_str("#FFFF99"));
            ("Y", "Tyr", Gadfly.@colorant_str("#FF66CC"));
            ("V", "Val", Gadfly.@colorant_str("#993300"));
            ("U", "Sec", Gadfly.@colorant_str("#00FFFF"));
            ("O", "Pyl", Gadfly.@colorant_str("#FFFF00"));
            ("B", "Asx", Gadfly.@colorant_str("#E4B9A3"));
            ("Z", "Glx", Gadfly.@colorant_str("#8F413C"));
            ("J", "Xle", Gadfly.@colorant_str("#825E12"));
            ("X", "Xaa", Gadfly.@colorant_str("#000000"));
            ("-",   "-", Gadfly.@colorant_str("#E0E0E0"));
            ("*",   "*", Gadfly.@colorant_str("#000000"));
            ( "",    "", Gadfly.@colorant_str("#FFFFFF"))]

  coltable = Dict(map(a -> (a[1], a[2]), colset))

  label = map(a -> a[2], colset)
  color = map(a -> a[3], colset)

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

  z = fill("", state.k * M)

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
        val = z[sub2ind((state.k, M), i, j)]

        (imax, jmax) = expand(i, j, state.k, M, val, z, processed)

        push!(xmin, j - 1)
        push!(xmax, jmax)
        push!(ymin, v[i])
        push!(ymax, v[imax + 1])
        push!(col, val)

        processed[i:imax, j:jmax] = true
        i = imax + 1
      else
        while i <= state.k && processed[i, j]
          i += 1
        end
      end
    end

    j += 1
  end

  idx = findin(label, col)

  majorfontsize = computefontsize(width, height)
  minorfontsize = 0.8 * majorfontsize

  keytitlefontsize = majorfontsize
  keylabelfontsize = 0.8 * keytitlefontsize
  padding = 1.0 * majorfontsize

  Gadfly.plot(Gadfly.layer(x=zeros(Float64, state.k - 1),
                           y=v[2:state.k],
                           xend=fill(float(M), state.k - 1),
                           yend=v[2:state.k],
                           Gadfly.Geom.segment),
              Gadfly.layer(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                           color=col, Gadfly.Geom.rect),
              Gadfly.Coord.cartesian(yflip=true,
                                     xmin=0, xmax=M,
                                     ymin=0, ymax=1),
              Gadfly.Scale.x_continuous(minvalue=1, maxvalue=M),
              Gadfly.Scale.y_continuous(minvalue=0, maxvalue=1,
                                        labels=l ->
                string(clusterlabel[yax .== l][1])),
              Gadfly.Scale.color_discrete_manual(color[idx]...,
                                                 levels=label[idx]),
              Gadfly.Theme(default_color=Gadfly.@colorant_str("black"),
                           background_color=Gadfly.@colorant_str("white"),
                           panel_fill=Gadfly.@colorant_str("white"),
                           grid_line_width=0Measures.mm,
                           line_width=linesep ? 0.3Measures.mm : 0Measures.mm,
                           lowlight_color=c->c,
                           lowlight_opacity=0.0,
                           major_label_font_size=majorfontsize,
                           minor_label_font_size=minorfontsize,
                           key_title_font_size=keytitlefontsize,
                           key_label_font_size=keylabelfontsize,
                           plot_padding=padding),
              Gadfly.Guide.yticks(ticks=yax),
              Gadfly.Guide.xlabel("Site"),
              Gadfly.Guide.ylabel("Cluster"),
              Gadfly.Guide.title(""),
              Gadfly.Guide.colorkey("Amino Acid"))
end
