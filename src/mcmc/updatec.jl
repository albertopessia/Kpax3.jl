# This file is part of Kpax3. License is MIT.
function updateC!(priorC::PriorColPartition,
                  mcmcobj::AminoAcidMCMC)
  mcmcobj.logpC = rpostpartitioncols!(mcmcobj.C, mcmcobj.cl, mcmcobj.v,
                                      mcmcobj.n1s, priorC)

  mcmcobj.loglik = 0.0

  lidx = 0
  for b in 1:size(mcmcobj.C, 2)
    for g in mcmcobj.cl
      lidx = mcmcobj.C[g, b] + 4 * (b - 1)
      mcmcobj.loglik += logmarglik(mcmcobj.n1s[g, b], mcmcobj.v[g],
                                   priorC.A[lidx], priorC.B[lidx])
    end
  end

  nothing
end
