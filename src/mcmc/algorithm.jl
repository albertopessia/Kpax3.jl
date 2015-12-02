# This file is part of Kpax3. License is MIT.

function kpax3mcmc(data::Array{UInt8, 2},
                   priorR::PriorRowPartition,
                   priorC::PriorColPartition,
                   settings::KSettings,
                   mcmcobj::AminoAcidMCMC)
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

        split!(ij, neighbours, S, data, priorR, priorC,
               settings, support, mcmcobj)
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

        merge!(ij, neighbours, S, data, priorR, priorC,
               settings, support, mcmcobj)
      end
    elseif opburnin[t] == 2

    elseif opburnin[t] == 3

    end

    if verbose && (t % verbosestep == 0)
      println("Burnin: iteration ", t, " done.")
    end
  end

end
