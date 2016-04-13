# This file is part of Kpax3. License is MIT.

function loglikmerge!(hi::Int,
                      hj::Int,
                      ni::Vector{Float64},
                      vi::Int,
                      priorC::AminoAcidPriorCol,
                      support::KSupport,
                      state::AminoAcidState)
  support.loglik = 0.0
  lidx = 0
  g = 0

  for b in 1:support.m
    for l in 1:(support.k - 1)
      g = support.cl[l]
      lidx = support.C[l, b] + 4 * (b - 1)
      support.loglik += logmarglik(state.n1s[g, b], state.v[g],
                                   priorC.A[lidx], priorC.B[lidx])
    end

    lidx = support.C[support.k, b] + 4 * (b - 1)
    support.loglik += logmarglik(ni[b], vi, priorC.A[lidx], priorC.B[lidx])
  end

  nothing
end
