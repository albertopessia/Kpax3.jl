# This file is part of Kpax3. License is MIT.

function splitloglik!(hi::Int,
                      priorC::AminoAcidPriorCol,
                      support::KSupport,
                      mcmcobj::AminoAcidMCMC)
  lidx = 0
  loglik = 0.0
  idx = 1
  for g in mcmcobj.cl
    if g != hi
      for b in 1:size(support.C, 2)
        lidx = support.C[idx, b] + 4 * (b - 1)
        loglik += logmarglik(mcmcobj.n1s[g, b], mcmcobj.v[g], priorC.A[lidx],
                             priorC.B[lidx])
      end
    else
      for b in 1:size(support.C, 2)
        lidx = support.C[idx, b] + 4 * (b - 1)
        loglik += logmarglik(support.ni[b], support.vi, priorC.A[lidx],
                             priorC.B[lidx])
      end
    end

    idx += 1
  end

  for b in 1:size(support.C, 2)
    lidx = support.C[idx, b] + 4 * (b - 1)
    loglik += logmarglik(support.nj[b], support.vj, priorC.A[lidx],
                         priorC.B[lidx])
  end

  loglik
end
