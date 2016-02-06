# This file is part of Kpax3. License is MIT.

function logliksplit!(hi::Int,
                      priorC::AminoAcidPriorCol,
                      support::KSupport,
                      mcmcobj::AminoAcidMCMC)
  support.loglik = 0.0
  l = 1
  lidx = 0

  for b in 1:size(support.C, 2)
    l = 1

    for g in mcmcobj.cl
      lidx = support.C[l, b] + 4 * (b - 1)

      if g != hi
        support.loglik += logmarglik(mcmcobj.n1s[g, b], mcmcobj.v[g],
                                     priorC.A[lidx], priorC.B[lidx])
      else
        support.loglik += logmarglik(support.ni[b], support.vi, priorC.A[lidx],
                                     priorC.B[lidx])
      end

      l += 1
    end

    lidx = support.C[l, b] + 4 * (b - 1)
    support.loglik += logmarglik(support.nj[b], support.vj, priorC.A[lidx],
                                 priorC.B[lidx])
  end

  nothing
end
