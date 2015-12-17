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
  k = length(mcmcobj.cl) + 1

  logω = [0.0; 0.0; log(k - 1.0) - log(k); -log(k)]

  initsupport!(ij, S, k, data, logω, priorC, support)

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
    for b in 1:size(data, 1)
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

  logprC, logpocC = simcsplit!(k, hi, logω, priorC, support, mcmcobj)

  distwm = Distributions.Beta(settings.parawm + support.vi,
                              settings.parawm + support.vj)

  logsrR = logpriorrowsplitratio(k, support.vi, support.vj, priorR)

  loglik = logliksplit!(hi, priorC, support, mcmcobj)

  ratio = exp(logsrR +
              logprC - mcmcobj.logprC +
              loglik - mcmcobj.loglik +
              Distributions.logpdf(distwm, ws) + mcmcobj.logpocC -
              Distributions.logpdf(settings.distws, ws) - lq - logpocC)

  if ratio >= 1 || (ratio > 0 &&
                    Distributions.rand(Distributions.Bernoulli(ratio)) == 1)
    hj = findfirst(!mcmcobj.filledcluster)

    if hj > 0
      idx = 0
      for g in mcmcobj.cl
        mcmcobj.C[g, 1] = support.C[idx += 1, 1]

        if g == hi
          mcmcobj.v[g] = support.vi
          mcmcobj.n1s[g, 1] = support.ni[1]
          mcmcobj.unit[g] = copy(support.ui[1:support.vi])
        end
      end

      mcmcobj.C[hj, 1] = support.C[idx += 1, 1]

      mcmcobj.filledcluster[hj] = true

      mcmcobj.v[hj] = support.vj
      mcmcobj.n1s[hj, 1] = support.nj[1]
      mcmcobj.unit[hj] = copy(support.uj[1:support.vj])

      for b in 2:size(data, 1)
        idx = 0
        for g in mcmcobj.cl
          mcmcobj.C[g, b] = support.C[idx += 1, b]
          mcmcobj.n1s[g, b] = g != hi ? mcmcobj.n1s[g, b] : support.ni[b]
        end
        mcmcobj.C[hj, b] = support.C[idx += 1, b]
        mcmcobj.n1s[hj, b] = support.nj[b]
      end
    else
      hj = k

      # reallocate memory
      len = min(settings.maxclust,
                size(data, 2) - length(mcmcobj.filledcluster)) - 1

      C = zeros(UInt8, k + len, size(data, 1))

      filledcluster = falses(k + len)
      v = zeros(Int, k + len)
      n1s = zeros(Float64, k + len, size(data, 1))
      unit = Vector{Int}[zeros(Int, settings.maxunit) for g in 1:(k + len)]

      idx = 0
      for g in mcmcobj.cl
        C[g, 1] = support.C[idx += 1, 1]

        filledcluster[g] = true

        if g != hi
          v[g] = mcmcobj.v[g]
          n1s[g, 1] = mcmcobj.n1s[g, 1]
          unit[g] = copy(mcmcobj.unit[g])
        else
          v[g] = support.vi
          n1s[g, 1] = support.ni[1]
          unit[g] = copy(support.ui[1:support.vi])
        end
      end

      C[k, 1] = support.C[idx += 1, 1]

      filledcluster[k] = true

      v[k] = support.vj
      n1s[k, 1] = support.nj[1]
      unit[k] = copy(support.uj[1:support.vj])

      for b in 2:size(data, 1)
        idx = 0
        for g in mcmcobj.cl
          C[g, b] = support.C[idx += 1, b]
          n1s[g, b] = g != hi ? mcmcobj.n1s[g, b] : support.ni[b]
        end
        C[k, b] = support.C[idx += 1, b]
        n1s[k, b] = support.nj[b]
      end

      mcmcobj.C = C

      mcmcobj.filledcluster = filledcluster
      mcmcobj.cl = find(mcmcobj.filledcluster)

      mcmcobj.v = v
      mcmcobj.n1s = n1s
      mcmcobj.unit = unit
    end

    # move units to their new cluster
    mcmcobj.R[support.uj[1:support.vj]] = hj

    mcmcobj.logprR += logsrR
    mcmcobj.logprC = logprC
    mcmcobj.loglik = loglik

    mcmcobj.logpocC = logpocC

    priorC.ω = [1.0, 1.0, (k - 1.0) / k, 1.0 / k]
    priorC.logω = logω
  end

  nothing
end
