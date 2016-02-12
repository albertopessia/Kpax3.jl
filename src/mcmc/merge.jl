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
  k = length(mcmcobj.cl) - 1

  initsupport!(ij, S, k, data, priorC, support)

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
  ni = zeros(Float64, size(data, 1))
  for b in 1:size(data, 1)
    ni[b] = mcmcobj.n1s[hi, b] + mcmcobj.n1s[hj, b]
  end

  # allocate the neighbours of i and j
  for l in 1:S
    u = neighbours[l]
    lcp[1] = lcp[2] = 0.0

    # compute p(x_{u} | x_{hi,1:(u-1)}) and p(x_{u} | x_{hj,1:(u-1)})
    for b in 1:size(data, 1)
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
    performmerge!(hi, hj, ni, vi, priorC, support, mcmcobj)
  end

  nothing
end

function performmerge!(hi::Int,
                       hj::Int,
                       ni::Vector{Float64},
                       vi::Int,
                       priorC::PriorColPartition,
                       support::KSupport,
                       mcmcobj::AminoAcidMCMC)
  mcmcobj.R[mcmcobj.unit[hj][1:mcmcobj.v[hj]]] = hi

  mcmcobj.filledcluster[hj] = false
  mcmcobj.cl = find(mcmcobj.filledcluster)

  idx = 0
  for g in mcmcobj.cl
    mcmcobj.C[g, 1] = support.C[idx += 1, 1]

    if g == hi
      mcmcobj.unit[g] = [mcmcobj.unit[hi][1:mcmcobj.v[hi]];
                         mcmcobj.unit[hj][1:mcmcobj.v[hj]]]
      mcmcobj.n1s[g, 1] = ni[1]
      mcmcobj.v[g] = vi
    end
  end

  for b in 2:size(mcmcobj.C, 2)
    idx = 0
    for g in mcmcobj.cl
      mcmcobj.C[g, b] = support.C[idx += 1, b]

      if g == hi
        mcmcobj.n1s[g, b] = ni[b]
      end
    end
  end

  mcmcobj.logpR += support.lograR
  copy!(mcmcobj.logpC, support.logpC)
  mcmcobj.loglik = support.loglik

  copy!(priorC.logω, support.logω)

  nothing
end
