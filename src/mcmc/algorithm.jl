# This file is part of Kpax3. License is MIT.

function randompermute!(x::Vector{Int},
                        S::Int)
  for i = S:-1:2
    j = rand(1:i)
    x[i], x[j] = x[j], x[i]
  end

  nothing
end

function splitmerge!(ij::Vector{Int},
                     neighbours::Vector{Int},
                     data::Matrix{UInt8},
                     priorR::PriorRowPartition,
                     priorC::PriorColPartition,
                     settings::KSettings,
                     support::KSupport,
                     state::AminoAcidState)
  # cluster founders (units i and j)
  StatsBase.sample!(1:support.n, ij, replace=false, ordered=false)

  # clusters of i and j respectively
  gi = state.R[ij[1]]
  gj = state.R[ij[2]]

  # total number of neighbours
  S = 0

  # generic unit
  u = 0

  if gi == gj
    for l in 1:state.v[gi]
      u = state.unit[gi][l]
      if (u != ij[1]) && (u != ij[2])
        neighbours[S += 1] = u
      end
    end

    randompermute!(neighbours, S)

    split!(ij, neighbours, S, data, priorR, priorC, settings, support, state)
  else
    for l in 1:state.v[gi]
      u = state.unit[gi][l]
      if u != ij[1]
        neighbours[S += 1] = u
      end
    end

    for l in 1:state.v[gj]
      u = state.unit[gj][l]
      if u != ij[2]
        neighbours[S += 1] = u
      end
    end

    randompermute!(neighbours, S)

    merge!(ij, neighbours, S, data, priorR, priorC, settings, support, state)
  end

  nothing
end

function kpax3mcmc!(data::Matrix{UInt8},
                    priorR::PriorRowPartition,
                    priorC::PriorColPartition,
                    settings::KSettings,
                    support::KSupport,
                    state::AminoAcidState)
  fp = open(settings.ofile, "w")

  # indices of units i and j
  ij = zeros(Int, 2)

  # neighbour indices
  neighbours = zeros(Int, support.n)

  try
    write(fp, settings.α)
    write(fp, settings.θ)
    write(fp, settings.γ[1])
    write(fp, settings.γ[2])
    write(fp, settings.γ[3])
    write(fp, settings.r)

    write(fp, support.n)
    write(fp, support.m)

    if settings.burnin > 0
      if settings.verbose
        println("Starting burnin phase...")
      end

      # sample which operators we are going to use
      operator = StatsBase.sample([0x01; 0x02; 0x03], settings.op,
                                  settings.burnin)

      for t in 1:settings.burnin
        if operator[t] == 0x01
          splitmerge!(ij, neighbours, data, priorR, priorC, settings, support,
                      state)
        elseif operator[t] == 0x02
          biased_random_walk!(data, priorR, priorC, settings, support, state)
        elseif operator[t] == 0x03
          updateC!(priorC, state)
        end

        if settings.verbose && (t % settings.verbosestep == 0)
          println("Burnin: step ", t, " done.")
        end
      end

      if settings.verbose
        println("Burnin phase completed.")
      end
    end

    if settings.verbose
      println("Starting collecting samples...")
    end

    operator = StatsBase.sample([0x01; 0x02; 0x03], settings.op, settings.T)

    for t in 1:settings.T
      if operator[t] == 0x01
        splitmerge!(ij, neighbours, data, priorR, priorC, settings, support,
                    state)
      elseif operator[t] == 0x02
        biased_random_walk!(data, priorR, priorC, settings, support, state)
      elseif operator[t] == 0x03
        updateC!(priorC, state)
      end

      if t % settings.tstep == 0
        savestate!(fp, state)
      end

      if settings.verbose && (t % settings.verbosestep == 0)
        println("Step ", t, " done.")
      end
    end

    if settings.verbose
      println("Markov Chain simulation complete.")
    end
  finally
    close(fp)
  end

  nothing
end

function kpax3mcmc(settings::KSettings)
  R = normalizepartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(k, settings.maxclust))

  support = KSupport(size(x.data, 1), size(x.data, 2), settings.maxclust,
                     settings.maxunit)

  state = AminoAcidState(x.data, R, priorR, priorC, settings)

  kpax3mcmc!(x.data, priorR, priorC, settings, support, state)

  nothing
end

function kpax3mcmc(x::AminoAcidData,
                   partition,
                   settings::KSettings)
  R = normalizepartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(k, settings.maxclust))

  support = KSupport(size(x.data, 1), size(x.data, 2), settings.maxclust,
                     settings.maxunit)

  state = AminoAcidState(x.data, R, priorR, priorC, settings)

  kpax3mcmc!(x.data, priorR, priorC, settings, support, state)

  nothing
end
