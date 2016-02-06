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

  filledcluster::BitArray{1}
  cl::Vector{Int}

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

  filledcluster = falses(settings.maxclust)
  filledcluster[unique(R)] = true

  cl = find(filledcluster)

  v = zeros(Int, settings.maxclust)
  n1s = zeros(Float64, settings.maxclust, m)
  unit = Vector{Int}[zeros(Int, settings.maxunit) for g in 1:settings.maxclust]

  for a in 1:n
    g = R[a]

    if v[g] == length(unit[g])
      tmp = zeros(Int, min(v[g] + settings.maxunit, n))
      unit[g] = copy!(tmp, unit[g])
    end

    v[g] += 1
    unit[g][v[g]] = a

    for b in 1:m
      n1s[g, b] += Float64(data[b, a])
    end
  end

  logpR = logdPriorRow(n, length(cl), v, priorR)
  logpC = rpostpartitioncols!(C, cl, v, n1s, priorC)
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
    for g in cl
      lidx = C[g, b] + 4 * (b - 1)
      loglik += logmarglik(n1s[g, b], v[g], priorC.A[lidx], priorC.B[lidx])
    end
  end

  AminoAcidMCMC(R, C, filledcluster, cl, v, n1s, unit, logpR, logpC, loglik)
end
