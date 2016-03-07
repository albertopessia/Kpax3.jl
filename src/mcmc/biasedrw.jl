# This file is part of Kpax3. License is MIT.

function biased_random_walk!(data::Matrix{UInt8},
                             priorR::PriorRowPartition,
                             priorC::PriorColPartition,
                             settings::KSettings,
                             support::KSupport,
                             mcmcobj::AminoAcidMCMC)
  i = StatsBase.sample(1:size(data, 2))

  # "new" cluster will get the same label if mcmcobj.v[hi] == 1
  hi = mcmcobj.R[i]
  hj = hi

  k = length(mcmcobj.cl)
  v = StatsBase.sample(1:k)

  if v == k
    # move i to a new cluster
    if mcmcobj.v[hi] > 1
      k += 1
      hj = findfirst(!mcmcobj.filledcluster)
      support.lograR = logratiopriorrowbrwsplit(k, mcmcobj.v[hi], priorR)
    end
  else
    hj = v < hi ? mcmcobj.cl[v] : mcmcobj.cl[v + 1]

    if mcmcobj.v[hi] > 1
      support.lograR = logratiopriorrowbrwmove(mcmcobj.v[hi], mcmcobj.v[hj],
                                               priorR)
    else
      k -= 1
      support.lograR = logratiopriorrowbrwmerge(k, mcmcobj.v[hj], priorR)
    end
  end

  initsupportbrw!(k, i, mcmcobj.v[hi], data, support)

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
  if hi == hj
    performbrwupdate!(support, mcmcobj)
  else
    # remove i from the list of units of cluster hi
    idx = 0
    for j in 1:mcmcobj.v[hi]
      if mcmcobj.unit[hi][j] != i
        support.ui[idx += 1] = mcmcobj.unit[hi][j]
      end
    end

    if hj > 0
      if mcmcobj.v[hi] > 1
        if mcmcobj.filledcluster[hj]
          performbrwmove!(i, hi, hj, support, mcmcobj)
        else
          performbrwsplit!(i, hi, hj, priorC, support, mcmcobj)
        end
      else
        performbrwmerge!(i, hi, hj, priorC, support, mcmcobj)
      end
    else
      performbrwsplitallocate!(i, hi, k, priorC, settings, support, mcmcobj)
    end
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
                          support::KSupport,
                          mcmcobj::AminoAcidMCMC)
  mcmcobj.R[i] = hj

  mcmcobj.filledcluster[hi] = false
  mcmcobj.cl = find(mcmcobj.filledcluster)

  mcmcobj.v[hj] += 1

  push!(mcmcobj.unit[hj], i)

  for b in 1:size(mcmcobj.C, 2)
    idx = 0
    for g in mcmcobj.cl
      mcmcobj.C[g, b] = support.C[idx += 1, b]
    end

    mcmcobj.n1s[hj, b] += support.ni[b]
  end

  copy!(priorC.logω, support.logω)

  nothing
end

function performbrwmove!(i::Int,
                         hi::Int,
                         hj::Int,
                         support::KSupport,
                         mcmcobj::AminoAcidMCMC)
  mcmcobj.R[i] = hj

  mcmcobj.v[hi] -= 1
  mcmcobj.v[hj] += 1

  mcmcobj.unit[hi] = copy(support.ui[1:(support.vi - 1)])
  mcmcobj.unit[hj] = vcat(mcmcobj.unit[hj], i)

  for b in 1:size(mcmcobj.C, 2)
    idx = 0
    for g in mcmcobj.cl
      mcmcobj.C[g, b] = support.C[idx += 1, b]
    end

    mcmcobj.n1s[hi, b] -= support.ni[b]
    mcmcobj.n1s[hj, b] += support.ni[b]
  end

  nothing
end

function performbrwsplit!(i::Int,
                          hi::Int,
                          hj::Int,
                          priorC::PriorColPartition,
                          support::KSupport,
                          mcmcobj::AminoAcidMCMC)
  mcmcobj.R[i] = hj

  for b in 1:size(mcmcobj.C, 2)
    idx = 0
    for g in mcmcobj.cl
      mcmcobj.C[g, b] = support.C[idx += 1, b]
    end
    mcmcobj.C[hj, b] = support.C[idx += 1, b]

    mcmcobj.n1s[hi, b] -= support.ni[b]
    mcmcobj.n1s[hj, b] = support.ni[b]
  end

  mcmcobj.filledcluster[hj] = true
  mcmcobj.cl = find(mcmcobj.filledcluster)

  mcmcobj.v[hi] -= 1
  mcmcobj.v[hj] = 1

  mcmcobj.unit[hi] = copy(support.ui[1:(support.vi - 1)])
  mcmcobj.unit[hj] = [i]

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
  len = min(length(mcmcobj.R), k + settings.maxclust - 1)

  C = zeros(UInt8, len, size(mcmcobj.C, 2))

  filledcluster = falses(len)
  v = zeros(Int, len)
  n1s = zeros(Float64, len, size(mcmcobj.C, 2))
  unit = Vector{Int}[zeros(Int, 0) for g in 1:len]

  idx = 0
  for g in mcmcobj.cl
    C[g, 1] = support.C[idx += 1, 1]

    v[g] = mcmcobj.v[g]
    n1s[g, 1] = mcmcobj.n1s[g, 1]
    unit[g] = copy(mcmcobj.unit[g])

    filledcluster[g] = true
  end

  C[k, 1] = support.C[idx += 1, 1]

  filledcluster[k] = true

  v[hi] -= 1
  n1s[hi, 1] -= support.ni[1]
  unit[hi] = copy(support.ui[1:(support.vi - 1)])

  v[k] = 1
  n1s[k, 1] = support.ni[1]
  unit[k] = [i]

  for b in 2:size(mcmcobj.C, 2)
    idx = 0
    for g in mcmcobj.cl
      C[g, b] = support.C[idx += 1, b]
      n1s[g, b] = mcmcobj.n1s[g, b]
    end
    C[k, b] = support.C[idx += 1, b]

    n1s[hi, b] -= support.ni[b]
    n1s[k, b] = support.ni[b]
  end

  mcmcobj.R[i] = k

  mcmcobj.C = C

  mcmcobj.filledcluster = filledcluster
  mcmcobj.cl = find(mcmcobj.filledcluster)

  mcmcobj.v = v
  mcmcobj.n1s = n1s
  mcmcobj.unit = unit

  copy!(priorC.logω, support.logω)

  nothing
end

function performbrwupdate!(support::KSupport,
                           mcmcobj::AminoAcidMCMC)
  for b in 1:size(mcmcobj.C, 2)
    idx = 0
    for g in mcmcobj.cl
      mcmcobj.C[g, b] = support.C[idx += 1, b]
    end
  end

  nothing
end
