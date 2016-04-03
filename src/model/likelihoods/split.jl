# This file is part of Kpax3. License is MIT.

function logliksplit!(hi::Int,
                      priorC::AminoAcidPriorCol,
                      support::KSupport,
                      mcmcobj::AminoAcidMCMC)
  support.loglik = 0.0
  lidx = 0
  g = 0

  for b in 1:support.m
    for l in 1:(support.k - 2)
      g = support.cl[l]
      lidx = support.C[l, b] + 4 * (b - 1)
      support.loglik += logmarglik(mcmcobj.n1s[g, b], mcmcobj.v[g],
                                   priorC.A[lidx], priorC.B[lidx])
    end

    lidx = support.C[support.k - 1, b] + 4 * (b - 1)
    support.loglik += logmarglik(support.ni[b], support.vi, priorC.A[lidx],
                                 priorC.B[lidx])

    lidx = support.C[support.k, b] + 4 * (b - 1)
    support.loglik += logmarglik(support.nj[b], support.vj, priorC.A[lidx],
                                 priorC.B[lidx])
  end

  nothing
end
