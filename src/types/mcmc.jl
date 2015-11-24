# This file is part of Kpax3. License is MIT.

"""
# Control MCMC run

## Description

## Fields

* `R` Int vector representing row clusters
* `C` UInt8 matrix representing column clusters
* `T` Length of the Markov Chain (total number of iterations)
* `burnin` Length of the burnin period
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

  T::Int
  burnin::Int
  op::WeightVec
  t::Int
  outfile::AbstractString

  k::Int
  emptyclusters::BitArray{1}
  v::Array{Float64, 1}
  n1s::Array{Float64, 2}
end

function AminoAcidMCMC(data::Array{UInt8, 2},
                       R::Array{Int, 1},
                       priorC::AminoAcidPriorCol,
                       T::Int,
                       maxclust::Int,
                       burnin::Int,
                       op::Array{Float64, 1},
                       t::Int,
                       outfile::AbstractString)
  m, n = size(data)

  k = maximum(R)

  emptyclusters = trues(maxclust)
  emptyclusters[unique(R)] = false

  v = zeros(Float64, maxclust)
  n1s = zeros(Float64, m, maxclust)

  for i in 1:n
    cl = R[i]
    v[cl] += 1.0

    for j in 1:m
      n1s[j, cl] += data[j, i]
    end
  end

  C = zeros(UInt8, m, maxclust)
  rcolpartition!(priorC, C, v, n1s, emptyclusters)

  AminoAcidMCMC(R, C, T, burnin, WeightVec(op), t, outfile, k, emptyclusters, v, n1s)
end
