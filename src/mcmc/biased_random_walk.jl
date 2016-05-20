# This file is part of Kpax3. License is MIT.

function biased_random_walk!(data::Matrix{UInt8},
                             priorR::PriorRowPartition,
                             priorC::PriorColPartition,
                             settings::KSettings,
                             support::KSupport,
                             state::AminoAcidState)
  k = state.k

  i = sample(1:support.n)
  hi = state.R[i]

  if state.v[hi] > 1
    if k > 1
      v = sample(1:k)

      if v == k
        # move i to a new cluster
        k += 1
        hj = findfirst(state.emptycluster)
        support.lograR = logratiopriorrowbrwsplit(k, state.v[hi], priorR)
      else
        hj = state.cl[v] < hi ? state.cl[v] : state.cl[v + 1]
        support.lograR = logratiopriorrowbrwmove(state.v[hi], state.v[hj],
                                                 priorR)
      end
    else
      # move i to a new cluster
      k = 2
      hj = findfirst(state.emptycluster)
      support.lograR = logratiopriorrowbrwsplit(k, state.v[hi], priorR)
    end
  else
    # move i to another cluster
    k -= 1
    v = sample(1:k)
    hj = state.cl[v] < hi ? state.cl[v] : state.cl[v + 1]
    support.lograR = logratiopriorrowbrwmerge(k, state.v[hj], priorR)
  end

  initsupportbrw!(k, i, state.v[hi], data, priorC, settings, support)

  simcbrw!(k, hi, hj, priorC, support, state)

  loglikbrw!(k, hi, hj, priorC, support, state)

  ratio = exp(support.lograR +
              support.logpC[1] - state.logpC[1] +
              support.loglik - state.loglik +
              state.logpC[2] - support.logpC[2])

  if ratio >= 1 || ((ratio > 0) && (rand() <= ratio))
    performbrw!(i, hi, hj, k, priorC, settings, support, state)
  end

  nothing
end

function performbrw!(i::Int,
                     hi::Int,
                     hj::Int,
                     k::Int,
                     priorC::PriorColPartition,
                     settings::KSettings,
                     support::KSupport,
                     state::AminoAcidState)
  if hj == 0
    resizestate!(state, k, settings)
    hj = k
  end

  # remove i from the list of units of cluster hi
  idx = 0
  for j in 1:state.v[hi]
    if state.unit[hi][j] != i
      support.ui[idx += 1] = state.unit[hi][j]
    end
  end

  if state.v[hi] > 1
    if state.emptycluster[hj]
      performbrwsplit!(i, hi, hj, priorC, settings, support, state)
    else
      performbrwmove!(i, hi, hj, settings, support, state)
    end
  else
    performbrwmerge!(i, hi, hj, priorC, settings, support, state)
  end

  state.logpR += support.lograR
  copy!(state.logpC, support.logpC)
  state.loglik = support.loglik
  state.logpp = state.logpR + state.logpC[1] + state.loglik

  nothing
end

function performbrwmerge!(i::Int,
                          hi::Int,
                          hj::Int,
                          priorC::PriorColPartition,
                          settings::KSettings,
                          support::KSupport,
                          state::AminoAcidState)
  state.R[i] = hj

  state.emptycluster[hi] = true

  k = 0
  for a in 1:length(state.emptycluster)
    if !state.emptycluster[a]
      state.cl[k += 1] = a
    end
  end

  state.k = k

  state.v[hj] += 1

  if length(state.unit[hj]) < state.v[hj]
    tmp = zeros(Int, min(support.n, state.v[hj] + settings.maxunit - 1))
    state.unit[hj] = copy!(tmp, state.unit[hj])
  end

  state.unit[hj][state.v[hj]] = i

  for b in 1:support.m
    for l in 1:(support.k - 1)
      state.C[support.cl[l], b] = support.C[l, b]
    end
    state.C[hj, b] = support.C[support.k, b]

    state.n1s[hj, b] += support.ni[b]
  end

  nothing
end

function performbrwmove!(i::Int,
                         hi::Int,
                         hj::Int,
                         settings::KSettings,
                         support::KSupport,
                         state::AminoAcidState)
  state.R[i] = hj

  state.v[hi] -= 1
  state.v[hj] += 1

  copy!(state.unit[hi], 1, support.ui, 1, support.vi - 1)

  if length(state.unit[hj]) < state.v[hj]
    tmp = zeros(Int, min(support.n, state.v[hj] + settings.maxunit - 1))
    state.unit[hj] = copy!(tmp, state.unit[hj])
  end

  state.unit[hj][state.v[hj]] = i

  for b in 1:support.m
    for l in 1:(support.k - 2)
      state.C[support.cl[l], b] = support.C[l, b]
    end
    state.C[hi, b] = support.C[support.k - 1, b]
    state.C[hj, b] = support.C[support.k, b]

    state.n1s[hi, b] -= support.ni[b]
    state.n1s[hj, b] += support.ni[b]
  end

  nothing
end

function performbrwsplit!(i::Int,
                          hi::Int,
                          hj::Int,
                          priorC::PriorColPartition,
                          settings::KSettings,
                          support::KSupport,
                          state::AminoAcidState)
  state.R[i] = hj

  for b in 1:support.m
    for l in 1:(support.k - 2)
      state.C[support.cl[l], b] = support.C[l, b]
    end
    state.C[hi, b] = support.C[support.k - 1, b]
    state.C[hj, b] = support.C[support.k, b]

    state.n1s[hi, b] -= support.ni[b]
    state.n1s[hj, b] = support.ni[b]
  end

  state.emptycluster[hj] = false

  k = 0
  for a in 1:length(state.emptycluster)
    if !state.emptycluster[a]
      state.cl[k += 1] = a
    end
  end

  state.k = k

  state.v[hi] -= 1
  state.v[hj] = 1

  copy!(state.unit[hi], 1, support.ui, 1, support.vi - 1)

  if length(state.unit[hj]) < state.v[hj]
    tmp = zeros(Int, min(support.n, state.v[hj] + settings.maxunit - 1))
    state.unit[hj] = copy!(tmp, state.unit[hj])
  end

  state.unit[hj][state.v[hj]] = i

  nothing
end
