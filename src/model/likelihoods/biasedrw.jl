# This file is part of Kpax3. License is MIT.

function loglikbrw!(k::Int,
                    hi::Int,
                    hj::Int,
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
        if g != hj
          # this cluster is not being affected by the move
          support.loglik += logmarglik(mcmcobj.n1s[g, b], mcmcobj.v[g],
                                       priorC.A[lidx], priorC.B[lidx])
        else
          # we are moving i into an existing cluster
          support.loglik += logmarglik(mcmcobj.n1s[g, b] + support.ni[b],
                                       mcmcobj.v[g] + 1, priorC.A[lidx],
                                       priorC.B[lidx])
        end

        l += 1
      else
        if g == hj
          # move singleton i into the same cluster
          support.loglik += logmarglik(support.ni[b], 1, priorC.A[lidx],
                                       priorC.B[lidx])
          l += 1
        elseif mcmcobj.v[g] > 1
          # remove i from this cluster
          support.loglik += logmarglik(mcmcobj.n1s[g, b] - support.ni[b],
                                       mcmcobj.v[g] - 1, priorC.A[lidx],
                                       priorC.B[lidx])
          l += 1
        end
      end
    end

    if length(mcmcobj.cl) < k
      # move i into a new cluster of its own
      lidx = support.C[l, b] + 4 * (b - 1)
      support.loglik += logmarglik(support.ni[b], 1, priorC.A[lidx],
                                   priorC.B[lidx])
    end
  end

  nothing
end
