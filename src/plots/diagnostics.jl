# This file is part of Kpax3. License is MIT.

function plottrace(entropy::Vector{Float64};
                   maxlag::Int=200,
                   M::Int=20000,
                   main::String="",
                   width::Real=640.0,
                   height::Real=360.0)
  nsim = length(entropy)

  # plot the last M points
  if M < nsim
    xtr = max(1, nsim - M):nsim
    ytr = entropy[xtr]
  else
    xtr = 1:nsim
    ytr = copy(entropy)
  end

  title = if main != ""
    ac = StatsBase.autocov(entropy, 0:maxlag)

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

  Plots.plot(xtr, ytr,
             seriestype=:path, linecolor=:black,
             # plot options
             background_color=:white, html_output_format=:svg,
             size=(width, height), window_title="",
             # subplot options
             label="", legend=:none, title=title,
             # x axis
             xlabel="Iteration", xlims=(minimum(xtr), maximum(xtr)),
             # y axis
             ylabel="Entropy",
             ylims=(max(0.0, minimum(ytr) - 0.05), maximum(ytr) + 0.05))
end

function plotdensity(entropy::Vector{Float64};
                     maxlag::Int=200,
                     main::String="",
                     width::Real=640.0,
                     height::Real=360.0)
  nsim = length(entropy)

  title = if main != ""
    ac = StatsBase.autocov(entropy, 0:maxlag)

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

  StatPlots.density(entropy,
                    linecolor=:black,
                    # plot options
                    background_color=:white, html_output_format=:svg,
                    size=(width, height), window_title="",
                    # subplot options
                    label="", legend=:none, title=title,
                    # x axis
                    xlabel="Entropy",
                    xlims=(max(0.0, minimum(entropy) - 0.05),
                           maximum(entropy) + 0.05),
                    # y axis
                    ylabel="Density")
end

function plotjump(avgd::Vector{Float64};
                  main::String="",
                  width::Real=640.0,
                  height::Real=360.0)
  maxlag = length(avgd)

  # Exponential smoothing for estimating the long term mean
  asy = avgd[1]
  for t in 2:maxlag
    asy = 0.6 * asy + 0.4 * avgd[t]
  end

  p = Plots.plot(1:maxlag, avgd,
                 seriestype=:scatter, markercolor=:black,
                 # plot options
                 background_color=:white, html_output_format=:svg,
                 size=(width, height), window_title="",
                 # subplot options
                 label="", legend=:none, title=main,
                 # x axis
                 xlabel="Lag", xlims=(0, maxlag + 1),
                 # y axis
                 ylabel="Average distance")

  Plots.plot!(p, [asy], seriestype=:hline, linestyle=:dash, linecolor=:black)

  p
end

function plotdgn(entropy::Vector{Float64},
                 avgd::Vector{Float64};
                 maxlag::Int=200,
                 M::Int=20000,
                 main::String="",
                 width::Real=640.0,
                 height::Real=360.0)
  p1 = plottrace(entropy, maxlag=maxlag, M=M, main=main,
                 width=width, height=height)
  p2 = plotdensity(entropy, maxlag=maxlag, main="", width=width,
                   height=height)
  p3 = plotjump(avgd, main="", width=width, height=height)

  Plots.plot(p1, p2, p3,
             layout=Plots.@layout([a; b; c]),
             background_color=:white, html_output_format=:svg,
             size=(10.0 + width, 10.0 + 3 * height), window_title="",
             label="", legend=:none, title="")
end
