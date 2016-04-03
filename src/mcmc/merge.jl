# This file is part of Kpax3. License is MIT.

function merge!(ij::Vector{Int},
                neighbours::Vector{Int},
                S::Int,
                data::Matrix{UInt8},
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                settings::KSettings,
                support::KSupport,
                mcmcobj::AminoAcidMCMC)
  # number of clusters after the merge
  k = mcmcobj.k - 1

  initsupportsplitmerge!(ij, S, k, data, priorC, settings, support)

  hi = mcmcobj.R[ij[1]]
  hj = mcmcobj.R[ij[2]]

  distwm = Distributions.Beta(settings.parawm + mcmcobj.v[hi],
                              settings.parawm + mcmcobj.v[hj])

  # sample a new proportion for cluster 'hi'
  w = Distributions.rand(distwm)

  # logarithm of the product of sequential probabilities
  lq = 0.0

  # temporary / support variables
  u = 0
  lcp = zeros(Float64, 2)
  z = 0.0
  p = 0.0

  vi = mcmcobj.v[hi] + mcmcobj.v[hj]
  ni = zeros(Float64, support.m)
  for b in 1:support.m
    ni[b] = mcmcobj.n1s[hi, b] + mcmcobj.n1s[hj, b]
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

    if mcmcobj.R[u] == hi
      updateclusteri!(u, data, support)
      lq += z
    else
      updateclusterj!(u, data, support)
      lq += log1p(-p)
    end
  end

  simcmerge!(k, hi, hj, vi, ni, priorC, support, mcmcobj)

  support.lograR = logratiopriorrowmerge(k, support.vi, support.vj, priorR)

  loglikmerge!(hi, hj, ni, vi, priorC, support, mcmcobj)

  ratio = exp(support.lograR +
              support.logpC[1] - mcmcobj.logpC[1] +
              support.loglik - mcmcobj.loglik +
              mcmcobj.logpC[2] - support.logpC[2] +
              Distributions.logpdf(settings.distws, w) -
              Distributions.logpdf(distwm, w) +
              lq)

  if ratio >= 1 || ((ratio > 0) && (rand() <= ratio))
    performmerge!(hi, hj, ni, vi, priorC, settings, support, mcmcobj)
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
                       mcmcobj::AminoAcidMCMC)
  for j in 1:mcmcobj.v[hj]
    mcmcobj.R[mcmcobj.unit[hj][j]] = hi
  end

  mcmcobj.emptycluster[hj] = true

  k = 0
  for a in 1:length(mcmcobj.emptycluster)
    if !mcmcobj.emptycluster[a]
      mcmcobj.cl[k += 1] = a
    end
  end

  mcmcobj.k = k

  for b in 1:support.m
    for l in 1:(support.k - 1)
      mcmcobj.C[support.cl[l], b] = support.C[l, b]
    end
    mcmcobj.C[hi, b] = support.C[support.k, b]

    mcmcobj.n1s[hi, b] = ni[b]
  end

  if length(mcmcobj.unit[hi]) < vi
    tmp = zeros(Int, min(support.n, vi + settings.maxunit - 1))
    copy!(tmp, 1, mcmcobj.unit[hi], 1, mcmcobj.v[hi])
    copy!(tmp, mcmcobj.v[hi] + 1, mcmcobj.unit[hj], 1, mcmcobj.v[hj])
    mcmcobj.unit[hi] = tmp
  else
    copy!(mcmcobj.unit[hi], mcmcobj.v[hi] + 1, mcmcobj.unit[hj], 1,
          mcmcobj.v[hj])
  end

  mcmcobj.v[hi] = vi

  mcmcobj.logpR += support.lograR
  copy!(mcmcobj.logpC, support.logpC)
  mcmcobj.loglik = support.loglik

  copy!(priorC.logω, support.logω)

  nothing
end
