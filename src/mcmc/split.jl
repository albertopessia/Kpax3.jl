# This file is part of Kpax3. License is MIT.

function split!(ij::Vector{Int},
                neighbours::Vector{Int},
                S::Int,
                data::Matrix{UInt8},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport,
                mcmcobj::AminoAcidMCMC)
  # number of clusters after the split
  k = mcmcobj.k + 1

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

  hi = mcmcobj.R[ij[1]]

  simcsplit!(k, hi, priorC, support, mcmcobj)

  support.lograR = logratiopriorrowsplit(k, support.vi, support.vj, priorR)

  logliksplit!(hi, priorC, support, mcmcobj)

  distwm = Distributions.Beta(settings.parawm + support.vi,
                              settings.parawm + support.vj)

  ratio = exp(support.lograR +
              support.logpC[1] - mcmcobj.logpC[1] +
              support.loglik - mcmcobj.loglik +
              mcmcobj.logpC[2] - support.logpC[2] +
              Distributions.logpdf(distwm, w) -
              Distributions.logpdf(settings.distws, w) -
              lq)

  if ratio >= 1 || ((ratio > 0) && (rand() <= ratio))
    performsplit!(hi, k, priorC, settings, support, mcmcobj)
  end

  nothing
end

function performsplit!(hi::Int,
                       k::Int,
                       priorC::PriorColPartition,
                       settings::KSettings,
                       support::KSupport,
                       mcmcobj::AminoAcidMCMC)
  hj = findfirst(mcmcobj.emptycluster)

  if hj > 0
    for b in 1:support.m
      for l in 1:(support.k - 2)
        mcmcobj.C[support.cl[l], b] = support.C[l, b]
      end
      mcmcobj.C[hi, b] = support.C[support.k - 1, b]
      mcmcobj.C[hj, b] = support.C[support.k, b]

      mcmcobj.n1s[hi, b] = support.ni[b]
      mcmcobj.n1s[hj, b] = support.nj[b]
    end

    mcmcobj.emptycluster[hj] = false

    h = 0
    for a in 1:length(mcmcobj.emptycluster)
      if !mcmcobj.emptycluster[a]
        mcmcobj.cl[h += 1] = a
      end
    end

    mcmcobj.k = h

    mcmcobj.v[hi] = support.vi
    mcmcobj.unit[hi] = copy!(mcmcobj.unit[hi], 1, support.ui, 1, support.vi)

    mcmcobj.v[hj] = support.vj

    if length(mcmcobj.unit[hj]) < mcmcobj.v[hj]
      tmp = zeros(Int, min(support.n, mcmcobj.v[hj] + settings.maxunit - 1))
      mcmcobj.unit[hj] = copy!(tmp, 1, support.uj, 1, support.vj)
    else
      mcmcobj.unit[hj] = copy!(mcmcobj.unit[hj], 1, support.uj, 1, support.vj)
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
    unit = Vector{Int}[zeros(Int, 0) for g in 1:len]

    g = 0
    for l in 1:(support.k - 2)
      g = support.cl[l]
      C[g, 1] = support.C[l, 1]
      v[g] = mcmcobj.v[g]
      n1s[g, 1] = mcmcobj.n1s[g, 1]
      unit[g] = mcmcobj.unit[g]
      emptycluster[g] = false
    end

    C[hi, 1] = support.C[support.k - 1, 1]
    v[hi] = support.vi
    n1s[hi, 1] = support.ni[1]
    unit[hi] = copy!(mcmcobj.unit[hi], 1, support.ui, 1, support.vi)
    emptycluster[hi] = false

    C[k, 1] = support.C[support.k, 1]
    v[k] = support.vj
    n1s[k, 1] = support.nj[1]
    unit[k] = zeros(Int, min(support.n, v[k] + settings.maxunit - 1))
    copy!(unit[k], 1, support.uj, 1, support.vj)
    emptycluster[k] = false

    for b in 2:support.m
      for l in 1:(support.k - 2)
        g = support.cl[l]
        C[g, b] = support.C[l, b]
        n1s[g, b] = mcmcobj.n1s[g, b]
      end
      C[hi, b] = support.C[support.k - 1, b]
      n1s[hi, b] = support.ni[b]
      C[k, b] = support.C[support.k, b]
      n1s[k, b] = support.nj[b]
    end

    mcmcobj.C = C

    mcmcobj.emptycluster = emptycluster

    h = 0
    for a in 1:length(mcmcobj.emptycluster)
      if !mcmcobj.emptycluster[a]
        cl[h += 1] = a
      end
    end

    mcmcobj.cl = cl
    mcmcobj.k = h

    mcmcobj.v = v
    mcmcobj.n1s = n1s
    mcmcobj.unit = unit
  end

  # move units to their new cluster
  for j in 1:support.vj
    mcmcobj.R[support.uj[j]] = hj
  end

  mcmcobj.logpR += support.lograR
  copy!(mcmcobj.logpC, support.logpC)
  mcmcobj.loglik = support.loglik

  copy!(priorC.logω, support.logω)

  nothing
end
