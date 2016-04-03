# This file is part of Kpax3. License is MIT.
function updateC!(priorC::PriorColPartition,
                  mcmcobj::AminoAcidMCMC)
  mcmcobj.logpC = rpostpartitioncols!(mcmcobj.C, mcmcobj.cl, mcmcobj.k,
                                      mcmcobj.v, mcmcobj.n1s, priorC)

  mcmcobj.loglik = 0.0

  lidx = 0
  for b in 1:size(mcmcobj.C, 2)
    for l in 1:mcmcobj.k
      lidx = mcmcobj.C[mcmcobj.cl[l], b] + 4 * (b - 1)
      mcmcobj.loglik += logmarglik(mcmcobj.n1s[mcmcobj.cl[l], b],
                                   mcmcobj.v[mcmcobj.cl[l]], priorC.A[lidx],
                                   priorC.B[lidx])
    end
  end

  nothing
end
