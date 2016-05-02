# This file is part of Kpax3. License is MIT.

function merge!(ij::Vector{Int},
                neighbours::Vector{Int},
                S::Int,
                data::Matrix{UInt8},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport,
                state::AminoAcidState)
  # number of clusters after the merge
  k = state.k - 1

  initsupportsplitmerge!(ij, S, k, data, priorC, settings, support)

  hi = state.R[ij[1]]
  hj = state.R[ij[2]]

  distwm = Distributions.Beta(settings.parawm + state.v[hi],
                              settings.parawm + state.v[hj])

  # sample a new proportion for cluster 'hi'
  w = Distributions.rand(distwm)

  # logarithm of the product of sequential probabilities
  lq = 0.0

  # temporary / support variables
  u = 0
  lcp = zeros(Float64, 2)
  z = 0.0
  p = 0.0

  vi = state.v[hi] + state.v[hj]
  ni = zeros(Float64, support.m)
  for b in 1:support.m
    ni[b] = state.n1s[hi, b] + state.n1s[hj, b]
  end

  # allocate the neighbours of i and j
  for l in 1:S
    u = neighbours[l]
    lcp[1] = lcp[2] = 0.0

    # compute p(x_{u} | x_{hi,1:(u-1)}) and p(x_{u} | x_{hj,1:(u-1)})
    for b in 1:support.m
      lcp[1] += computeclusteriseqprobs!(data[b, u], b, priorC, support)
      lcp[2] += computeclusterjseqprobs!(data[b, u], b, priorC, support)
    end

    # (w * p1) / (w * p1 + (1 - w) * p2) = 1 / (1 + ((1 - w) / w) * (p2 / p1))
    # => e^(-log(1 + e^(log(1 - w) - log(w) + log(p2) - log(p1))))
    z = -log1p(exp(log(1 - w) - log(w) + lcp[2] - lcp[1]))
    p = exp(z)

    if state.R[u] == hi
      updateclusteri!(u, data, support)
      lq += z
    else
      updateclusterj!(u, data, support)
      lq += log1p(-p)
    end
  end

  simcmerge!(k, hi, hj, vi, ni, priorC, support, state)

  support.lograR = logratiopriorrowmerge(k, support.vi, support.vj, priorR)

  loglikmerge!(hi, hj, ni, vi, priorC, support, state)

  ratio = exp(support.lograR +
              support.logpC[1] - state.logpC[1] +
              support.loglik - state.loglik +
              state.logpC[2] - support.logpC[2] +
              Distributions.logpdf(settings.distws, w) -
              Distributions.logpdf(distwm, w) +
              lq)

  if ratio >= 1 || ((ratio > 0) && (rand() <= ratio))
    performmerge!(hi, hj, ni, vi, priorC, settings, support, state)
  end

  nothing
end

function performmerge!(hi::Int,
                       hj::Int,
                       ni::Vector{Float64},
                       vi::Int,
                       priorC::PriorColPartition,
                       settings::KSettings,
                       support::KSupport,
                       state::AminoAcidState)
  for j in 1:state.v[hj]
    state.R[state.unit[hj][j]] = hi
  end

  state.emptycluster[hj] = true

  k = 0
  for a in 1:length(state.emptycluster)
    if !state.emptycluster[a]
      state.cl[k += 1] = a
    end
  end

  state.k = k

  for b in 1:support.m
    for l in 1:(support.k - 1)
      state.C[support.cl[l], b] = support.C[l, b]
    end
    state.C[hi, b] = support.C[support.k, b]

    state.n1s[hi, b] = ni[b]
  end

  if length(state.unit[hi]) < vi
    tmp = zeros(Int, min(support.n, vi + settings.maxunit - 1))
    copy!(tmp, 1, state.unit[hi], 1, state.v[hi])
    copy!(tmp, state.v[hi] + 1, state.unit[hj], 1, state.v[hj])
    state.unit[hi] = tmp
  else
    copy!(state.unit[hi], state.v[hi] + 1, state.unit[hj], 1,
          state.v[hj])
  end

  state.v[hi] = vi

  state.logpR += support.lograR
  copy!(state.logpC, support.logpC)
  state.loglik = support.loglik

  nothing
end
