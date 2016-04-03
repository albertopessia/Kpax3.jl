# This file is part of Kpax3. License is MIT.

function loglikbrw!(k::Int,
                    hi::Int,
                    hj::Int,
                    priorC::AminoAcidPriorCol,
                    support::KSupport,
                    mcmcobj::AminoAcidMCMC)
  support.loglik = 0.0
  lidx = 0
  g = 0

  if k >= mcmcobj.k
    h = support.k - 2
  else
    h = support.k - 1
  end

  for b in 1:support.m
    for l in 1:h
      g = support.cl[l]
      lidx = support.C[l, b] + 4 * (b - 1)
      support.loglik += logmarglik(mcmcobj.n1s[g, b], mcmcobj.v[g],
                                   priorC.A[lidx], priorC.B[lidx])
    end

    if k == mcmcobj.k
      lidx = support.C[support.k - 1, b] + 4 * (b - 1)
      support.loglik += logmarglik(mcmcobj.n1s[hi, b] - support.ni[b],
                                   mcmcobj.v[hi] - 1, priorC.A[lidx],
                                   priorC.B[lidx])
      lidx = support.C[support.k, b] + 4 * (b - 1)
      support.loglik += logmarglik(mcmcobj.n1s[hj, b] + support.ni[b],
                                   mcmcobj.v[hj] + 1, priorC.A[lidx],
                                   priorC.B[lidx])
    elseif k > mcmcobj.k
      lidx = support.C[support.k - 1, b] + 4 * (b - 1)
      support.loglik += logmarglik(mcmcobj.n1s[hi, b] - support.ni[b],
                                   mcmcobj.v[hi] - 1, priorC.A[lidx],
                                   priorC.B[lidx])
     lidx = support.C[support.k, b] + 4 * (b - 1)
     support.loglik += logmarglik(support.ni[b], 1, priorC.A[lidx],
                                  priorC.B[lidx])
    else
      lidx = support.C[support.k, b] + 4 * (b - 1)
      support.loglik += logmarglik(mcmcobj.n1s[hj, b] + support.ni[b],
                                   mcmcobj.v[hj] + 1, priorC.A[lidx],
                                   priorC.B[lidx])
    end
  end

  nothing
end
