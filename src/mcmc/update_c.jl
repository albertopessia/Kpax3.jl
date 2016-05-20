# This file is part of Kpax3. License is MIT.

function updateC!(priorC::PriorColPartition,
                  state::AminoAcidState)
  state.logpC = rpostpartitioncols!(state.C, state.cl, state.k, state.v,
                                    state.n1s, priorC)

  state.loglik = loglikelihood(state.C, state.cl, state.k, state.v, state.n1s,
                               priorC)

  state.logpp = state.logpR + state.logpC[1] + state.loglik

  nothing
end
