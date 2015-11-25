# This file is part of Kpax3. License is MIT.

function kpax3mcmc(mcmcobj::AminoAcidMCMC,
                   priorR::PriorRowPartition,
                   priorC::PriorColPartition,
                   data::Array{UInt8, 2},
                   verbose::Bool,
                   verbosestep::Int)
  m, n = size(data)

  opburnin = sample(1:3, mcmcobj.op, mcmcobj.burnin)
  opmcmc = sample(1:3, mcmcobj.op, mcmcobj.T)

  ij = zeros(Int, 2)
  neighbours = falses(n)

  if verbose && (mcmcobj.burnin > 0)
    println("Starting burnin phase...")
  end

  for t in 1:length(opburnin)
    if opburnin[t] == 1
      # cluster founders (units i and j)
      sample!(1:n, ij; replace=false, ordered=false)

      # neighbours of units i and j
      neighbours[:] = (mcmcobj.R .== mcmcobj.R[ij[1]]) |
                      (mcmcobj.R .== mcmcobj.R[ij[2]])
      neighbours[ij] = false

      S = find(neighbours)

      if mcmcobj.R[ij[1]] == mcmcobj.R[ij[2]]
        split!(mcmcobj, priorR, priorC, data, ij, S)
      else
        merge!(mcmcobj, priorR, priorC, data, ij, S)
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
