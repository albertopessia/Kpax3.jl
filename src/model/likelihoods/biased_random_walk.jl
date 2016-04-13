# This file is part of Kpax3. License is MIT.

function loglikbrw!(k::Int,
                    hi::Int,
                    hj::Int,
                    priorC::AminoAcidPriorCol,
                    support::KSupport,
                    state::AminoAcidState)
  support.loglik = 0.0
  lidx = 0
  g = 0

  if k >= state.k
    h = support.k - 2
  else
    h = support.k - 1
  end

  for b in 1:support.m
    for l in 1:h
      g = support.cl[l]
      lidx = support.C[l, b] + 4 * (b - 1)
      support.loglik += logmarglik(state.n1s[g, b], state.v[g],
                                   priorC.A[lidx], priorC.B[lidx])
    end

    if k == state.k
      lidx = support.C[support.k - 1, b] + 4 * (b - 1)
      support.loglik += logmarglik(state.n1s[hi, b] - support.ni[b],
                                   state.v[hi] - 1, priorC.A[lidx],
                                   priorC.B[lidx])
      lidx = support.C[support.k, b] + 4 * (b - 1)
      support.loglik += logmarglik(state.n1s[hj, b] + support.ni[b],
                                   state.v[hj] + 1, priorC.A[lidx],
                                   priorC.B[lidx])
    elseif k > state.k
      lidx = support.C[support.k - 1, b] + 4 * (b - 1)
      support.loglik += logmarglik(state.n1s[hi, b] - support.ni[b],
                                   state.v[hi] - 1, priorC.A[lidx],
                                   priorC.B[lidx])
     lidx = support.C[support.k, b] + 4 * (b - 1)
     support.loglik += logmarglik(support.ni[b], 1, priorC.A[lidx],
                                  priorC.B[lidx])
    else
      lidx = support.C[support.k, b] + 4 * (b - 1)
      support.loglik += logmarglik(state.n1s[hj, b] + support.ni[b],
                                   state.v[hj] + 1, priorC.A[lidx],
                                   priorC.B[lidx])
    end
  end

  nothing
end
