# This file is part of Kpax3. License is MIT.

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
type AminoAcidMCMC <: Kpax3MCMC
  R::Array{Int, 1}
  C::Array{UInt8, 2}
  k::Int
  emptycluster::BitArray{1}
  cluster::Array{Cluster, 1}
end

function AminoAcidMCMC(data::Array{UInt8, 2},
                       R::Array{Int, 1},
                       priorC::AminoAcidPriorCol,
                       settings::Kpax3Settings)
  m, n = size(data)

  k = maximum(R)

  C = zeros(UInt8, m, settings.maxclust)

  emptycluster = trues(settings.maxclust)
  emptycluster[unique(R)] = false

  cluster = [Cluster(0, zeros(Int, 1), zeros(Float64, 1), 0.0)
             for g in 1:settings.maxclust]

  for g in 1:k
    cluster[g].unit = zeros(Int, settings.maxunit)
    cluster[g].n1s = zeros(Float64, m)
  end

  for i in 1:n
    g = R[i]
    cluster[g].v += 1

    if cluster[g].v > length(cluster[g].unit)
      tmp = zeros(Int, min(2 * settings.maxunit, n))
      tmp[1:(cluster[g].v - 1)] = cluster[g].unit
      cluster[g].unit = tmp
    end

    cluster[g].unit[cluster[g].v] = i

    for j in 1:m
      cluster[g].n1s[j] += data[j, i]
    end
  end

  rcolpartition!(priorC, C, cluster, emptycluster)

  # If array A has dimension (d_{1}, ..., d_{l}, ..., d_{L}), to access
  # element A[i_{1}, ..., i_{l}, ..., i_{L}] it is possible to use the
  # following linear index
  # linearidx = i_{1} + d_{1} * (i_{2} - 1) + ... +
  #           + (d_{1} * ... * d_{l-1}) * (i_{l} - 1) + ... +
  #           + (d_{1} * ... * d_{L-1}) * (i_{L} - 1)
  #
  # A[i_{1}, ..., i_{l}, ..., i_{L}] == A[linearidx]
  for g in 1:k
    linearidx = [(i + m * (C[i, g] - 1))::Int for i in 1:m]
    cluster[g].ll = sum(logmarglik(cluster[g].n1s, cluster[g].v,
                                   priorC.A[linearidx], priorC.B[linearidx]))
  end

  AminoAcidMCMC(R, C, k, emptycluster, cluster)
end
