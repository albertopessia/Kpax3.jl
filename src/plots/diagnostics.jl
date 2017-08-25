# This file is part of Kpax3. License is MIT.

function plottrace(entropy::Vector{Float64};
                   ac::Vector{Float64}=zeros(Float64, 0),
                   maxlag::Int=200,
                   M::Int=20000,
                   main::String="",
                   width::Real=183.0,
                   height::Real=92.0)
  nsim = length(entropy)

  # plot the last M points
  if M < nsim
    xtr = max(1, nsim - M):nsim
    ytr = entropy[xtr]
  else
    xtr = 1:nsim
    ytr = StatsBase.deepcopy(entropy)
  end

  majorfontsize = computefontsize(width, height)
  minorfontsize = 0.8 * majorfontsize

  keytitlefontsize = majorfontsize
  keylabelfontsize = 0.8 * keytitlefontsize
  linewidth = 0.05 * majorfontsize
  padding = 1.0 * majorfontsize
  gridwidth = 0.1 * majorfontsize

  title = if main != ""
    if length(ac) == 0
      ac = StatsBase.autocov(entropy, 0:maxlag)
    end

    variid = ac[1] / nsim
    mcvar = imsevar(ac, nsim)
    mcse = sqrt(mcvar)
    eff = min(nsim, ess(variid, mcvar, nsim))

    effstr = if eff <= 100000
      @sprintf("%d", eff)
    else
      @sprintf("%.3e", eff)
    end

    sestr = if 0.001 <= mcse <= 1000
      @sprintf("%.3f", mcse)
    else
      @sprintf("%.3e", mcse)
    end

    string(main, " (mcse = ", sestr, ", ESS = ", effstr, ")")
  else
    ""
  end

  Gadfly.plot(Gadfly.layer(x=xtr, y=ytr, Gadfly.Geom.line),
              Gadfly.Coord.cartesian(ymin=max(0.0, minimum(ytr) - 0.05)),
              Gadfly.Theme(default_color=Gadfly.@colorant_str("black"),
                           background_color=Gadfly.@colorant_str("white"),
                           panel_fill=Gadfly.@colorant_str("white"),
                           major_label_font_size=majorfontsize,
                           minor_label_font_size=minorfontsize,
                           key_title_font_size=keytitlefontsize,
                           key_label_font_size=keylabelfontsize,
                           plot_padding=padding,
                           line_width=linewidth,
                           grid_line_width=gridwidth),
              Gadfly.Guide.xlabel("Iteration"),
              Gadfly.Guide.ylabel("Entropy"),
              Gadfly.Guide.title(title))
end

function plotdensity(entropy::Vector{Float64};
                     ac::Vector{Float64}=zeros(Float64, 0),
                     maxlag::Int=200,
                     main::String="",
                     width::Real=183.0,
                     height::Real=92.0)
  nsim = length(entropy)

  majorfontsize = computefontsize(width, height)
  minorfontsize = 0.8 * majorfontsize

  keytitlefontsize = majorfontsize
  keylabelfontsize = 0.8 * keytitlefontsize
  linewidth = 0.05 * majorfontsize
  padding = 1.0 * majorfontsize
  gridwidth = 0.1 * majorfontsize

  title = if main != ""
    if length(ac) == 0
      ac = StatsBase.autocov(entropy, 0:maxlag)
    end

    variid = ac[1] / nsim
    mcvar = imsevar(ac, nsim)
    mcse = sqrt(mcvar)
    eff = min(nsim, ess(variid, mcvar, nsim))

    effstr = if eff <= 100000
      @sprintf("%d", eff)
    else
      @sprintf("%.3e", eff)
    end

    sestr = if 0.001 <= mcse <= 1000
      @sprintf("%.3f", mcse)
    else
      @sprintf("%.3e", mcse)
    end

    string(main, " (mcse = ", sestr, ", ESS = ", effstr, ")")
  else
    ""
  end

  Gadfly.plot(Gadfly.layer(x=entropy, Gadfly.Geom.density),
              Gadfly.Coord.cartesian(xmin=max(0.0, minimum(entropy) - 0.1)),
              Gadfly.Theme(default_color=Gadfly.@colorant_str("black"),
                           background_color=Gadfly.@colorant_str("white"),
                           panel_fill=Gadfly.@colorant_str("white"),
                           major_label_font_size=majorfontsize,
                           minor_label_font_size=minorfontsize,
                           key_title_font_size=keytitlefontsize,
                           key_label_font_size=keylabelfontsize,
                           plot_padding=padding,
                           line_width=linewidth,
                           grid_line_width=gridwidth),
              Gadfly.Guide.xlabel("Entropy"),
              Gadfly.Guide.ylabel("Density"),
              Gadfly.Guide.title(title))
end

function plotjump(avgd::Vector{Float64};
                  main::String="",
                  width::Real=183.0,
                  height::Real=92.0)
  maxlag = length(avgd)

  majorfontsize = computefontsize(width, height)
  minorfontsize = 0.8 * majorfontsize

  keytitlefontsize = majorfontsize
  keylabelfontsize = 0.8 * keytitlefontsize
  pointsize = 0.3 * majorfontsize
  linewidth = 0.05 * majorfontsize
  padding = 1.0 * majorfontsize
  gridwidth = 0.1 * majorfontsize

  asy = avgd[1]
  for t in 2:maxlag
    asy = 0.6 * asy + 0.4 * avgd[t]
  end

  xticks = Gadfly.optimize_ticks(1, maxlag, strict_span=true)[1]
  if xticks[1] == 0.0
    xticks[1] = 1.0
  end

  if xticks[end] > maxlag
    xticks[end] = maxlag
  elseif xticks[end] < maxlag
    push!(xticks, maxlag)
  end

  Gadfly.plot(Gadfly.layer(x=1:maxlag, y=avgd, Gadfly.Geom.point),
              Gadfly.layer(yintercept=[asy],
                           Gadfly.Geom.hline(size=linewidth)),
              Gadfly.Coord.cartesian(xmin=0, xmax=maxlag + 1),
              Gadfly.Theme(default_color=Gadfly.@colorant_str("black"),
                           background_color=Gadfly.@colorant_str("white"),
                           panel_fill=Gadfly.@colorant_str("white"),
                           major_label_font_size=majorfontsize,
                           minor_label_font_size=minorfontsize,
                           key_title_font_size=keytitlefontsize,
                           key_label_font_size=keylabelfontsize,
                           point_size=pointsize,
                           plot_padding=padding,
                           line_width=linewidth,
                           grid_line_width=gridwidth),
              Gadfly.Guide.xticks(ticks=xticks),
              Gadfly.Guide.xlabel("Lag"),
              Gadfly.Guide.ylabel("Average distance"),
              Gadfly.Guide.title(main))
end

function plotdgn(entropy::Vector{Float64},
                 avgd::Vector{Float64};
                 ac::Vector{Float64}=zeros(Float64, 0),
                 maxlag::Int=200,
                 M::Int=20000,
                 main::String="",
                 width::Real=183.0,
                 height::Real=92.0)
  p1 = plottrace(entropy, ac=ac, maxlag=maxlag, M=M, main=main, width=width,
                 height=height)
  p2 = plotdensity(entropy, ac=ac, maxlag=maxlag, main="", width=width,
                   height=height)
  p3 = plotjump(avgd, main="", width=width, height=height)

  Gadfly.vstack(p1, p2, p3)
end
