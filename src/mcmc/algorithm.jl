# This file is part of Kpax3. License is MIT.

function randompermute!(x::Array{Int, 1},
                        S::Int)
  for i = S:-1:2
    j = rand(1:i)
    x[i], x[j] = x[j], x[i]
  end

  nothing
end

function saveresults!(fp::IOStream,
                      mcmcobj::AminoAcidMCMC)
  write(fp, mcmcobj.k)
  write(fp, mcmcobj.R)
  write(fp, vec(mcmcobj.C[!mcmcobj.emptycluster, :]))

  nothing
end

function splitmerge!(ij::Array{Int, 1},
                     neighbours::Array{Int, 1},
                     data::Array{UInt8, 2},
                     priorR::PriorRowPartition,
                     priorC::PriorColPartition,
                     settings::KSettings,
                     support::KSupport,
                     mcmcobj::AminoAcidMCMC)
  # cluster founders (units i and j)
  StatsBase.sample!(1:size(data, 2), ij, replace=false, ordered=false)

  # clusters of i and j respectively
  gi = mcmcobj.R[ij[1]]
  gj = mcmcobj.R[ij[2]]

  # total number of neighbours
  S = zero(Int)

  if gi == gj
    for u in mcmcobj.cluster[gi].unit[1:mcmcobj.cluster[gi].v]
      if (u != ij[1]) && (u != ij[2])
        neighbours[S += 1] = u
      end
    end

    randompermute!(neighbours, S)

    split!(ij, neighbours, S, data, priorR, priorC, settings, support, mcmcobj)
  else
    for u in mcmcobj.cluster[gi].unit[1:mcmcobj.cluster[gi].v]
      if u != ij[1]
        neighbours[S += 1] = u
      end
    end

    for u in mcmcobj.cluster[gj].unit[1:mcmcobj.cluster[gj].v]
      if u != ij[2]
        neighbours[S += 1] = u
      end
    end

    randompermute!(neighbours, S)

    merge!(ij, neighbours, S, data, priorR, priorC, settings, support, mcmcobj)
  end

  nothing
end

function kpax3mcmc!(data::Array{UInt8, 2},
                    priorR::PriorRowPartition,
                    priorC::PriorColPartition,
                    settings::KSettings,
                    support::KSupport,
                    mcmcobj::AminoAcidMCMC)
  fp = open(settings.outfile, "w")

  # indices of units i and j
  ij = zeros(Int, 2)

  # neighbour indices
  neighbours = zeros(Int, size(data, 2))

  try
    write(fp, size(data, 2))
    write(fp, size(data, 1))

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
                      mcmcobj)
        elseif operator[t] == 0x02
          nothing
        elseif operator[t] == 0x03
          nothing
        end

        if settings.verbose && (t % settings.verbosestep == 0)
          println("Burnin: iteration ", t, " done.")
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
                    mcmcobj)
      elseif operator[t] == 0x02
        nothing
      elseif operator[t] == 0x03
        nothing
      end

      if t % settings.tstep == 0
        saveresults!(fp, mcmcobj)
      end

      if settings.verbose && (t % settings.verbosestep == 0)
        println("Iteration ", t, " done.")
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
