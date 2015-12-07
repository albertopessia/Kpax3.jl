# This file is part of Kpax3. License is MIT.

abstract KMCMC

"""
# Control MCMC run

## Description

## Fields

* `R` Int vector representing row clusters
* `C` UInt8 matrix representing column clusters
* `op` MCMC operators probabilities
* `t` Save the Markov Chain state `(R, C)` every `t` iterations
* `outfile` Path to the output file
* `k` Current number of clusters (number of distinct values in `R`)
* `emptyclusters` Bool vector representing which clusters are empty
* `v` Vector of cluster sizes
* `counts` Matrix containing the total number of ones observed in each cluster
"""
type AminoAcidMCMC <: KMCMC
  R::Array{Int, 1}
  C::Array{UInt8, 2}

  cluster::Array{KCluster, 1}
  emptycluster::BitArray{1}
  k::Int

  logprR::Float64
  logprC::Float64
  loglik::Float64

  logpocC::Float64
end

function AminoAcidMCMC(data::Array{UInt8, 2},
                       R::Array{Int, 1},
                       priorR::PriorRowPartition,
                       priorC::AminoAcidPriorCol,
                       settings::KSettings)
  m, n = size(data)
  k = maximum(R)

  C = zeros(UInt8, settings.maxclust, m)

  emptycluster = trues(settings.maxclust)
  emptycluster[unique(R)] = false

  cluster = [KCluster(0, zeros(Int, 1), zeros(Float64, 1))
             for g in 1:settings.maxclust]

  for g in 1:k
    cluster[g].unit = zeros(Int, settings.maxunit)
    cluster[g].n1s = zeros(Float64, m)
  end

  for i in 1:n
    g = R[i]
    cluster[g].v += 1

    if cluster[g].v > length(cluster[g].unit)
      tmp = zeros(Int, min(cluster[g].v - 1 + settings.maxunit, n))
      tmp[1:(cluster[g].v - 1)] = cluster[g].unit
      cluster[g].unit = tmp
    end

    cluster[g].unit[cluster[g].v] = i

    for j in 1:m
      cluster[g].n1s[j] += Float64(data[j, i])
    end
  end

  rpostpartitioncols!(C, cluster, emptycluster, priorC.logγ, priorC.logω,
                      priorC.A, priorC.B)

  logprR = logdPriorRow(n, k, [Int(cluster[g].v) for g in 1:k], priorR)
  logprC = logpriorC(C, emptycluster, priorC.logγ, priorC.logω)
  loglik = zero(Float64)

  for g in 1:k
    # If array A has dimension (d_{1}, ..., d_{l}, ..., d_{L}), to access
    # element A[i_{1}, ..., i_{l}, ..., i_{L}] it is possible to use the
    # following linear index
    # linearidx = i_{1} + d_{1} * (i_{2} - 1) + ... +
    #           + (d_{1} * ... * d_{l-1}) * (i_{l} - 1) + ... +
    #           + (d_{1} * ... * d_{L-1}) * (i_{L} - 1)
    #
    # A[i_{1}, ..., i_{l}, ..., i_{L}] == A[linearidx]
    linearidx = [(C[g, b] + 4 * (b - 1))::Int for b in 1:m]
    loglik += sum(logmarglik(cluster[g].n1s, cluster[g].v, priorC.A[linearidx],
                             priorC.B[linearidx]))
  end

  logpocC = logcondpostC(C, cluster, emptycluster, priorC.logγ, priorC.logω,
                         priorC.A, priorC.B)

  AminoAcidMCMC(R, C, cluster, emptycluster, k, logprR, logprC, loglik, logpocC)
end
