# This file is part of Kpax3. License is MIT.

function split!(ij::Vector{Int},
                neighbours::Vector{Int},
                S::Int,
                data::Matrix{UInt8},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport,
                state::AminoAcidState)
  # number of clusters after the split
  k = state.k + 1

  initsupportsplitmerge!(ij, S, k, data, priorC, settings, support)

  # sample a new proportion for cluster 'hi'
  w = Distributions.rand(settings.distws)

  # logarithm of the product of sequential probabilities
  lq = 0.0

  # temporary / support variables
  u = 0
  lcp = zeros(Float64, 2)
  z = 0.0
  p = 0.0

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

    if rand() <= p
      updateclusteri!(u, data, support)
      lq += z
    else
      updateclusterj!(u, data, support)
      lq += log1p(-p)
    end
  end

  hi = state.R[ij[1]]

  simcsplit!(k, hi, priorC, support, state)

  support.lograR = logratiopriorrowsplit(k, support.vi, support.vj, priorR)

  logliksplit!(hi, priorC, support, state)

  distwm = Distributions.Beta(settings.parawm + support.vi,
                              settings.parawm + support.vj)

  ratio = exp(support.lograR +
              support.logpC[1] - state.logpC[1] +
              support.loglik - state.loglik +
              state.logpC[2] - support.logpC[2] +
              Distributions.logpdf(distwm, w) -
              Distributions.logpdf(settings.distws, w) -
              lq)

  if ratio >= 1 || ((ratio > 0) && (rand() <= ratio))
    performsplit!(hi, k, priorC, settings, support, state)
  end

  nothing
end

function performsplit!(hi::Int,
                       k::Int,
                       priorC::PriorColPartition,
                       settings::KSettings,
                       support::KSupport,
                       state::AminoAcidState)
  hj = findfirst(state.emptycluster)

  if hj > 0
    for b in 1:support.m
      for l in 1:(support.k - 2)
        state.C[support.cl[l], b] = support.C[l, b]
      end
      state.C[hi, b] = support.C[support.k - 1, b]
      state.C[hj, b] = support.C[support.k, b]

      state.n1s[hi, b] = support.ni[b]
      state.n1s[hj, b] = support.nj[b]
    end

    state.emptycluster[hj] = false

    h = 0
    for a in 1:length(state.emptycluster)
      if !state.emptycluster[a]
        state.cl[h += 1] = a
      end
    end

    state.k = h

    state.v[hi] = support.vi
    state.unit[hi] = copy!(state.unit[hi], 1, support.ui, 1, support.vi)

    state.v[hj] = support.vj

    if length(state.unit[hj]) < state.v[hj]
      tmp = zeros(Int, min(support.n, state.v[hj] + settings.maxunit - 1))
      state.unit[hj] = copy!(tmp, 1, support.uj, 1, support.vj)
    else
      state.unit[hj] = copy!(state.unit[hj], 1, support.uj, 1, support.vj)
    end
  else
    hj = k

    # reallocate memory
    len = min(support.n, k + settings.maxclust - 1)

    C = zeros(UInt8, len, support.m)
    emptycluster = trues(len)
    cl = zeros(Int, len)
    v = zeros(Int, len)
    n1s = zeros(Float64, len, support.m)
    unit = Array{Vector{Int}}(len)

    # prevent losing pre-allocated vectors
    for l in 1:length(state.unit)
      unit[l] = state.unit[l]
    end

    for l in (length(state.unit) + 1):len
      unit[l] = zeros(Int, settings.maxunit)
    end

    g = 0
    for l in 1:(support.k - 2)
      g = support.cl[l]
      C[g, 1] = support.C[l, 1]
      v[g] = state.v[g]
      n1s[g, 1] = state.n1s[g, 1]
      emptycluster[g] = false
    end

    C[hi, 1] = support.C[support.k - 1, 1]
    v[hi] = support.vi
    n1s[hi, 1] = support.ni[1]
    copy!(state.unit[hi], 1, support.ui, 1, support.vi)
    emptycluster[hi] = false

    C[k, 1] = support.C[support.k, 1]
    v[k] = support.vj
    n1s[k, 1] = support.nj[1]

    if v[k] > length(unit[k])
      resize!(unit[k], min(support.n, v[k] + settings.maxunit - 1))
    end

    copy!(unit[k], 1, support.uj, 1, support.vj)
    emptycluster[k] = false

    for b in 2:support.m
      for l in 1:(support.k - 2)
        g = support.cl[l]
        C[g, b] = support.C[l, b]
        n1s[g, b] = state.n1s[g, b]
      end
      C[hi, b] = support.C[support.k - 1, b]
      n1s[hi, b] = support.ni[b]
      C[k, b] = support.C[support.k, b]
      n1s[k, b] = support.nj[b]
    end

    state.C = C

    state.emptycluster = emptycluster

    h = 0
    for a in 1:length(state.emptycluster)
      if !state.emptycluster[a]
        cl[h += 1] = a
      end
    end

    state.cl = cl
    state.k = h

    state.v = v
    state.n1s = n1s
    state.unit = unit
  end

  # move units to their new cluster
  for j in 1:support.vj
    state.R[support.uj[j]] = hj
  end

  state.logpR += support.lograR
  copy!(state.logpC, support.logpC)
  state.loglik = support.loglik

  nothing
end
