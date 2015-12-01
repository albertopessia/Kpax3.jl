# This file is part of Kpax3. License is MIT.

function kpax3mcmc(mcmcobj::AminoAcidMCMC,
                   priorR::PriorRowPartition,
                   priorC::PriorColPartition,
                   data::Array{UInt8, 2},
                   settings::Kpax3Settings)
  m, n = size(data)

  # sample before hand which operators we are going to use
  opburnin = StatsBase.sample(1:3, settings.op, settings.burnin)
  opmcmc = StatsBase.sample(1:3, settings.op, settings.T)

  # indices of units i and j
  ij = zeros(Int, 2)

  # neighbour indices
  neighbours = zeros(Int, n)

  # total number of neighbours
  S = 0

  # clusters of i and j respectively
  gi = 0
  gj = 0

  if verbose && (settings.burnin > 0)
    println("Starting burnin phase...")
  end

  for t in 1:burnin
    if opburnin[t] == 1
      # cluster founders (units i and j)
      StatsBase.sample!(1:n, ij, replace=false, ordered=false)

      # ij = [303, 347];
      gi = mcmcobj.R[ij[1]]
      gj = mcmcobj.R[ij[2]]

      S = 0

      if gi == gj
        for u in 1:mcmcobj.cluster[gi].v
          idx = mcmcobj.cluster[gi].unit[u]

          if (idx != ij[1]) && (idx != ij[2])
            neighbours[S += 1] = idx
          end
        end

        split!(mcmcobj, priorR, priorC, data, ij, neighbours, S, support)
      else
        for u in 1:mcmcobj.cluster[gi].v
          idx = mcmcobj.cluster[gi].unit[u]

          if idx != ij[1]
            neighbours[S += 1] = idx
          end
        end
        for u in 1:mcmcobj.cluster[gj].v
          idx = mcmcobj.cluster[gj].unit[u]

          if idx != ij[2]
            neighbours[S += 1] = idx
          end
        end

        merge!(mcmcobj, priorR, priorC, data, ij, neighbours, S, support)
      end
    elseif opburnin[t] == 2
      biasedrandomwalk!(mcmcobj, priorR, priorC, data)
    elseif opburnin[t] == 3
      updateC!(mcmcobj, priorR, priorC, data)
    end

    if verbose && (t % verbosestep == 0)
      println("Burnin: iteration ", t, " done.")
    end
  end

end
