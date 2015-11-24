# This file is part of Kpax3. License is MIT.

function kpax3mcmc(mcmcobj::AminoAcidMCMC,
                   priorR::EwensPitman,
                   priorC::AminoAcidPriorCol,
                   data::AminoAcidData,
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

end
