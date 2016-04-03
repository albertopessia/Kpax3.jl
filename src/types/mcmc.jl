# This file is part of Kpax3. License is MIT.

abstract KMCMC

"""
# Control MCMC run

## Description

## Fields

"""
type AminoAcidMCMC <: KMCMC
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

function AminoAcidMCMC(data::Matrix{UInt8},
                       R::Vector{Int},
                       priorR::PriorRowPartition,
                       priorC::AminoAcidPriorCol,
                       settings::KSettings)
  m, n = size(data)

  C = zeros(UInt8, settings.maxclust, m)

  emptycluster = trues(settings.maxclust)
  cl = zeros(Int, settings.maxclust)

  v = zeros(Int, settings.maxclust)
  n1s = zeros(Float64, settings.maxclust, m)
  unit = Vector{Int}[zeros(Int, settings.maxunit) for g in 1:settings.maxclust]

  g = 0
  for a in 1:n
    g = R[a]

    if emptycluster[g]
      emptycluster[g] = false
    end

    if v[g] == length(unit[g])
      tmp = zeros(Int, min(v[g] + settings.maxunit, n))
      unit[g] = copy!(tmp, unit[g])
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
  k = 0
  for a in 1:length(emptycluster)
    if !emptycluster[a]
      cl[k += 1] = a
    end
  end

  logpR = logdPriorRow(n, k, v, priorR)
  logpC = rpostpartitioncols!(C, cl, k, v, n1s, priorC)
  loglik = 0.0

  # If array A has dimension (d_{1}, ..., d_{l}, ..., d_{L}), to access
  # element A[i_{1}, ..., i_{l}, ..., i_{L}] it is possible to use the
  # following linear index
  # lidx = i_{1} + d_{1} * (i_{2} - 1) + ... +
  #      + (d_{1} * ... * d_{l-1}) * (i_{l} - 1) + ... +
  #      + (d_{1} * ... * d_{L-1}) * (i_{L} - 1)
  #
  # A[i_{1}, ..., i_{l}, ..., i_{L}] == A[lidx]
  lidx = 0
  for b in 1:m
    for l in 1:k
      lidx = C[cl[l], b] + 4 * (b - 1)
      loglik += logmarglik(n1s[cl[l], b], v[cl[l]], priorC.A[lidx],
                           priorC.B[lidx])
    end
  end

  AminoAcidMCMC(R, C, emptycluster, cl, k, v, n1s, unit, logpR, logpC, loglik)
end
