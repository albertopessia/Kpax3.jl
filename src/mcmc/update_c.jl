# This file is part of Kpax3. License is MIT.

function updateC!(priorC::PriorColPartition,
                  state::AminoAcidState)
  state.logpC = rpostpartitioncols!(state.C, state.cl, state.k, state.v,
                                    state.n1s, priorC)

  state.loglik = loglikelihood(state.C, state.cl, state.k, state.v, state.n1s,
                               priorC)

  nothing
end
