# This file is part of Kpax3. License is MIT.

abstract State

type AminoAcidState <: State
  R::Vector{Int}
  C::Matrix{UInt8}

  emptycluster::BitArray{1}
  cl::Vector{Int}
  k::Int

  v::Vector{Int}
  n1s::Matrix{Float64}
  unit::Vector{Vector{Int}}

  logpR::Float64
  logpC::Vector{Float64}
  loglik::Float64
end

function AminoAcidState(data::Matrix{UInt8},
                        R::Vector{Int},
                        priorR::PriorRowPartition,
                        priorC::AminoAcidPriorCol,
                        settings::KSettings)
  (m, n) = size(data)

  R = normalizepartition(R, n)
  k = maximum(R)

  maxclust = max(k, min(n, settings.maxclust))

  C = zeros(UInt8, maxclust, m)

  emptycluster = trues(maxclust)
  cl = zeros(Int, maxclust)

  v = zeros(Int, maxclust)
  n1s = zeros(Float64, maxclust, m)
  unit = Vector{Int}[zeros(Int, settings.maxunit) for g in 1:maxclust]

  g = 0
  for a in 1:n
    g = R[a]

    if emptycluster[g]
      emptycluster[g] = false
    end

    if v[g] == length(unit[g])
      resize!(unit[g], min(n, v[g] + settings.maxunit))
    end

    v[g] += 1
    unit[g][v[g]] = a

    for b in 1:m
      n1s[g, b] += float(data[b, a])
    end
  end

  # we do it here and not in the previous loop because we want them sorted
  # it is not really necessary but they will be sorted during the simulation
  # anyway (we will loop on emptycluster from now on)
  i = 1
  for a in 1:length(emptycluster)
    if !emptycluster[a]
      cl[i] = a
      i += 1
    end
  end

  logpR = logdPriorRow(n, k, v, priorR)

  logpC = zeros(Float64, 2)
  computelocalmode!(v, n1s, C, cl, k, logpC, priorC)

  loglik = loglikelihood(C, cl, k, v, n1s, priorC)

  AminoAcidState(R, C, emptycluster, cl, k, v, n1s, unit, logpR, logpC, loglik)
end

function copystate(x::AminoAcidState)
  R = copy(x.R)
  C = copy(x.C)

  emptycluster = copy(x.emptycluster)
  cl = copy(x.cl)
  k = x.k
  v = copy(x.v)
  n1s = copy(x.n1s)

  unit = Array{Vector{Int}}(length(x.unit))
  for l in 1:length(x.unit)
    unit[l] = copy(x.unit[l])
  end

  logpR = x.logpR
  logpC = copy(x.logpC)
  loglik = x.loglik

  AminoAcidState(R, C, emptycluster, cl, k, v, n1s, unit, logpR, logpC, loglik)
end

function copystate!(dest::AminoAcidState,
                    src::AminoAcidState)
  dest.k = src.k
  dest.logpR = src.logpR
  dest.logpC[1] = src.logpC[1]
  dest.logpC[2] = src.logpC[2]
  dest.loglik = src.loglik

  if length(dest.R) == length(src.R)
    copy!(dest.R, src.R)
  else
    dest.R = copy(src.R)
  end

  if (size(dest.C, 2) == size(src.C, 2)) && (size(dest.C, 1) >= src.cl[src.k])
    fill!(dest.emptycluster, true)

    g = 0
    for l in 1:src.k
      g = src.cl[l]

      dest.C[g, 1] = src.C[g, 1]
      dest.emptycluster[g] = false
      dest.cl[l] = src.cl[l]
      dest.v[g] = src.v[g]
      dest.n1s[g, 1] = src.n1s[g, 1]

      if length(dest.unit[g]) >= src.v[g]
        copy!(dest.unit[g], 1, src.unit[g], 1, src.v[g])
      else
        dest.unit[g] = copy(src.unit[g])
      end
    end

    for b in 2:size(dest.C, 2)
      for l in 1:src.k
        g = src.cl[l]

        dest.C[g, b] = src.C[g, b]
        dest.n1s[g, b] = src.n1s[g, b]
      end
    end
  else
    dest.C = copy(src.C)

    dest.emptycluster = copy(src.emptycluster)
    dest.cl = copy(src.cl)
    dest.v = copy(src.v)
    dest.n1s = copy(src.n1s)

    dest.unit = Array{Vector{Int}}(length(src.unit))
    for l in 1:length(src.unit)
      dest.unit[l] = copy(src.unit[l])
    end
  end

  nothing
end

function initializestate(data::Matrix{UInt8},
                         D::Matrix{Float64},
                         kset::UnitRange{Int},
                         priorR::PriorRowPartition,
                         priorC::PriorColPartition,
                         settings::KSettings)
  if length(priorC.logω) < kset[end]
    resizelogω!(priorC, kset[end])
  end

  n = size(data, 2)

  R = ones(Int, n)

  s = AminoAcidState(data, R, priorR, priorC, settings)
  slp = s.logpR + s.logpC[1] + s.loglik

  t1 = copystate(s)
  tlp1 = slp

  t2 = copystate(s)
  tlp2 = slp

  # TODO: remove the hack once kmedoids is fixed

  if settings.verbose
    @printf("Log-posterior (plus a constant) for one cluster: %.4f.\n", slp)
    @printf("Now scanning %d to %d clusters.\n", kset[1], kset[end])
  end

  niter = 0
  notupdated = true
  for k in kset
    if settings.verbose && (k % 10 == 0)
      @printf("Total number of clusters = %d.\n", k)
    end

    fill!(R, 0)
    notupdated = true
    while notupdated
      try
        copy!(R, kmedoids(D, k).assignments)
      catch
        fill!(R, 0)
      end

      if R[1] > 0
        notupdated = false
      end
    end

    updatestate!(t1, data, R, priorR, priorC, settings)
    tlp1 = t1.logpR + t1.logpC[1] + t1.loglik

    niter = 0
    while niter < 10
      notupdated = true
      while notupdated
        try
          copy!(R, kmedoids(D, k).assignments)
        catch
          fill!(R, 0)
        end

        if R[1] > 0
          notupdated = false
        end
      end

      updatestate!(t2, data, R, priorR, priorC, settings)
      tlp2 = t2.logpR + t2.logpC[1] + t2.loglik

      if tlp2 > tlp1
        copystate!(t1, t2)
        tlp1 = tlp2
      end

      niter += 1
    end

    if tlp1 > slp
      copystate!(s, t1)
      slp = tlp1

      if settings.verbose
        @printf("Found a better solution! ")
        @printf("Log-posterior (plus a constant) for %d clusters: %.4f.\n", k,
                slp)
      end
    end
  end

  (s, slp)
end

function updatestate!(state::AminoAcidState,
                      data::Matrix{UInt8},
                      R::Vector{Int},
                      priorR::PriorRowPartition,
                      priorC::AminoAcidPriorCol,
                      settings::KSettings)
  (m, n) = size(data)

  copy!(state.R, normalizepartition(R, n))
  state.k = maximum(R)

  if size(state.C, 1) < state.k
    # reallocate memory
    len = min(n, state.k + settings.maxclust - 1)

    state.C = zeros(UInt8, len, m)
    state.emptycluster = trues(len)
    state.cl = zeros(Int, len)
    state.v = zeros(Int, len)
    state.n1s = zeros(Float64, len, m)

    unit = Array{Vector{Int}}(len)
    for l in 1:length(state.unit)
      unit[l] = state.unit[l]
    end
    for l in (length(state.unit) + 1):len
      unit[l] = zeros(Int, settings.maxunit)
    end
    state.unit = unit
  else
    # TODO: no need to fill the whole array
    fill!(state.emptycluster, true)
    fill!(state.v, 0)
    fill!(state.n1s, 0.0)
  end

  g = 0
  for a in 1:n
    g = state.R[a]

    if state.emptycluster[g]
      state.emptycluster[g] = false
    end

    if state.v[g] == length(state.unit[g])
      resize!(state.unit[g], min(n, state.v[g] + settings.maxunit))
    end

    state.v[g] += 1
    state.unit[g][state.v[g]] = a

    for b in 1:m
      state.n1s[g, b] += float(data[b, a])
    end
  end

  i = 1
  for a in 1:length(state.emptycluster)
    if !state.emptycluster[a]
      state.cl[i] = a
      i += 1
    end
  end

  state.logpR = logdPriorRow(n, state.k, state.v, priorR)

  computelocalmode!(state.v, state.n1s, state.C, state.cl, state.k, state.logpC,
                    priorC)

  state.loglik = loglikelihood(state.C, state.cl, state.k, state.v, state.n1s,
                               priorC)

  nothing
end
