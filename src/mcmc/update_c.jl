# This file is part of Kpax3. License is MIT.
function updateC!(priorC::PriorColPartition,
                  state::AminoAcidState)
  state.logpC = rpostpartitioncols!(state.C, state.cl, state.k,
                                      state.v, state.n1s, priorC)

  state.loglik = 0.0

  lidx = 0
  for b in 1:size(state.C, 2)
    for l in 1:state.k
      lidx = state.C[state.cl[l], b] + 4 * (b - 1)
      state.loglik += logmarglik(state.n1s[state.cl[l], b],
                                   state.v[state.cl[l]], priorC.A[lidx],
                                   priorC.B[lidx])
    end
  end

  nothing
end
