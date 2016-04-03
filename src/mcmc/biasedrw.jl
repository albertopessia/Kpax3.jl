# This file is part of Kpax3. License is MIT.

function biased_random_walk!(data::Matrix{UInt8},
                             priorR::PriorRowPartition,
                             priorC::PriorColPartition,
                             settings::KSettings,
                             support::KSupport,
                             mcmcobj::AminoAcidMCMC)
  k = mcmcobj.k

  i = StatsBase.sample(1:support.n)
  hi = mcmcobj.R[i]

  if mcmcobj.v[hi] > 1
    if k > 1
      v = StatsBase.sample(1:k)

      if v == k
        # move i to a new cluster
        k += 1
        hj = findfirst(mcmcobj.emptycluster)
        support.lograR = logratiopriorrowbrwsplit(k, mcmcobj.v[hi], priorR)
      else
        hj = mcmcobj.cl[v] < hi ? mcmcobj.cl[v] : mcmcobj.cl[v + 1]
        support.lograR = logratiopriorrowbrwmove(mcmcobj.v[hi], mcmcobj.v[hj],
                                                 priorR)
      end
    else
      # move i to a new cluster
      k = 2
      hj = findfirst(mcmcobj.emptycluster)
      support.lograR = logratiopriorrowbrwsplit(k, mcmcobj.v[hi], priorR)
    end
  else
    # move i to another cluster
    k -= 1
    v = StatsBase.sample(1:k)
    hj = mcmcobj.cl[v] < hi ? mcmcobj.cl[v] : mcmcobj.cl[v + 1]
    support.lograR = logratiopriorrowbrwmerge(k, mcmcobj.v[hj], priorR)
  end

  initsupportbrw!(k, i, mcmcobj.v[hi], data, settings, support)

  simcbrw!(k, hi, hj, priorC, support, mcmcobj)

  loglikbrw!(k, hi, hj, priorC, support, mcmcobj)

  ratio = exp(support.lograR +
              support.logpC[1] - mcmcobj.logpC[1] +
              support.loglik - mcmcobj.loglik +
              mcmcobj.logpC[2] - support.logpC[2])

  if ratio >= 1 || ((ratio > 0) && (rand() <= ratio))
    performbrw!(i, hi, hj, k, priorC, settings, support, mcmcobj)
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
                     mcmcobj::AminoAcidMCMC)
  # remove i from the list of units of cluster hi
  idx = 0
  for j in 1:mcmcobj.v[hi]
    if mcmcobj.unit[hi][j] != i
      support.ui[idx += 1] = mcmcobj.unit[hi][j]
    end
  end

  if hj > 0
    if mcmcobj.v[hi] > 1
      if mcmcobj.emptycluster[hj]
        performbrwsplit!(i, hi, hj, priorC, settings, support, mcmcobj)
      else
        performbrwmove!(i, hi, hj, settings, support, mcmcobj)
      end
    else
      performbrwmerge!(i, hi, hj, priorC, settings, support, mcmcobj)
    end
  else
    performbrwsplitallocate!(i, hi, k, priorC, settings, support, mcmcobj)
  end

  mcmcobj.logpR += support.lograR
  copy!(mcmcobj.logpC, support.logpC)
  mcmcobj.loglik = support.loglik

  nothing
end

function performbrwmerge!(i::Int,
                          hi::Int,
                          hj::Int,
                          priorC::PriorColPartition,
                          settings::KSettings,
                          support::KSupport,
                          mcmcobj::AminoAcidMCMC)
  mcmcobj.R[i] = hj

  mcmcobj.emptycluster[hi] = true

  k = 0
  for a in 1:length(mcmcobj.emptycluster)
    if !mcmcobj.emptycluster[a]
      mcmcobj.cl[k += 1] = a
    end
  end

  mcmcobj.k = k

  mcmcobj.v[hj] += 1

  if length(mcmcobj.unit[hj]) < mcmcobj.v[hj]
    tmp = zeros(Int, min(support.n, mcmcobj.v[hj] + settings.maxunit - 1))
    mcmcobj.unit[hj] = copy!(tmp, mcmcobj.unit[hj])
  end

  mcmcobj.unit[hj][mcmcobj.v[hj]] = i

  for b in 1:support.m
    for l in 1:(support.k - 1)
      mcmcobj.C[support.cl[l], b] = support.C[l, b]
    end
    mcmcobj.C[hj, b] = support.C[support.k, b]

    mcmcobj.n1s[hj, b] += support.ni[b]
  end

  copy!(priorC.logω, support.logω)

  nothing
end

function performbrwmove!(i::Int,
                         hi::Int,
                         hj::Int,
                         settings::KSettings,
                         support::KSupport,
                         mcmcobj::AminoAcidMCMC)
  mcmcobj.R[i] = hj

  mcmcobj.v[hi] -= 1
  mcmcobj.v[hj] += 1

  copy!(mcmcobj.unit[hi], 1, support.ui, 1, support.vi - 1)

  if length(mcmcobj.unit[hj]) < mcmcobj.v[hj]
    tmp = zeros(Int, min(support.n, mcmcobj.v[hj] + settings.maxunit - 1))
    mcmcobj.unit[hj] = copy!(tmp, mcmcobj.unit[hj])
  end

  mcmcobj.unit[hj][mcmcobj.v[hj]] = i

  for b in 1:support.m
    for l in 1:(support.k - 2)
      mcmcobj.C[support.cl[l], b] = support.C[l, b]
    end
    mcmcobj.C[hi, b] = support.C[support.k - 1, b]
    mcmcobj.C[hj, b] = support.C[support.k, b]

    mcmcobj.n1s[hi, b] -= support.ni[b]
    mcmcobj.n1s[hj, b] += support.ni[b]
  end

  nothing
end

function performbrwsplit!(i::Int,
                          hi::Int,
                          hj::Int,
                          priorC::PriorColPartition,
                          settings::KSettings,
                          support::KSupport,
                          mcmcobj::AminoAcidMCMC)
  mcmcobj.R[i] = hj

  for b in 1:support.m
    for l in 1:(support.k - 2)
      mcmcobj.C[support.cl[l], b] = support.C[l, b]
    end
    mcmcobj.C[hi, b] = support.C[support.k - 1, b]
    mcmcobj.C[hj, b] = support.C[support.k, b]

    mcmcobj.n1s[hi, b] -= support.ni[b]
    mcmcobj.n1s[hj, b] = support.ni[b]
  end

  mcmcobj.emptycluster[hj] = false

  k = 0
  for a in 1:length(mcmcobj.emptycluster)
    if !mcmcobj.emptycluster[a]
      mcmcobj.cl[k += 1] = a
    end
  end

  mcmcobj.k = k

  mcmcobj.v[hi] -= 1
  mcmcobj.v[hj] = 1

  copy!(mcmcobj.unit[hi], 1, support.ui, 1, support.vi - 1)

  if length(mcmcobj.unit[hj]) < mcmcobj.v[hj]
    tmp = zeros(Int, min(support.n, mcmcobj.v[hj] + settings.maxunit - 1))
    mcmcobj.unit[hj] = copy!(tmp, mcmcobj.unit[hj])
  end

  mcmcobj.unit[hj][mcmcobj.v[hj]] = i

  copy!(priorC.logω, support.logω)

  nothing
end

function performbrwsplitallocate!(i::Int,
                                  hi::Int,
                                  k::Int,
                                  priorC::PriorColPartition,
                                  settings::KSettings,
                                  support::KSupport,
                                  mcmcobj::AminoAcidMCMC)
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
  v[hi] = mcmcobj.v[hi] - 1
  n1s[hi, 1] = mcmcobj.n1s[hi, 1] - support.ni[1]
  unit[hi] = copy!(mcmcobj.unit[hi], 1, support.ui, 1, support.vi - 1)
  emptycluster[hi] = false

  C[k, 1] = support.C[support.k, 1]
  v[k] = 1
  n1s[k, 1] = support.ni[1]
  unit[k] = zeros(Int, settings.maxunit)
  unit[k][1] = i
  emptycluster[k] = false

  for b in 2:support.m
    for l in 1:(support.k - 2)
      g = support.cl[l]
      C[g, b] = support.C[l, b]
      n1s[g, b] = mcmcobj.n1s[g, b]
    end
    C[hi, b] = support.C[support.k - 1, b]
    n1s[hi, b] = mcmcobj.n1s[hi, b] - support.ni[b]
    C[k, b] = support.C[support.k, b]
    n1s[k, b] = support.ni[b]
  end

  mcmcobj.R[i] = k

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

  copy!(priorC.logω, support.logω)

  nothing
end
