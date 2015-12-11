# This file is part of Kpax3. License is MIT.

"""
# User defined settings for a Kpax3 run

## Description

## Fields

* `outfile` Path to the output file
* `T` Length of the Markov Chain starting after the burnin period
* `burnin` Length of the burnin period
* `t` Save the Markov Chain state `(R, C)` every `t` iterations
* `op` WeightVec representing the probabilities of Markov Chain kernels
* `α` Ewens-Pitman formula 'discount' parameter
* `θ` Ewens-Pitman formula 'concentration' parameter
* `γ` Status probabilities
* `r` Scalar used to define Beta distribution parameters
* `maxunit` Maximum number of units per cluster for an initial memory allocation
* `maxclust` Maximum number of clusters for an initial memory allocation
* `verbose` If `true`, print status reports
* `verbosestep` Print a status report every `verbosestep` Markov Chain steps

"""
immutable KSettings
  outfile::AbstractString
  T::Int
  burnin::Int
  tstep::Int
  op::StatsBase.WeightVec
  α::Float64
  θ::Float64
  γ::Vector{Float64}
  r::Float64
  distws::Distributions.Beta
  parawm::Float64
  maxclust::Int
  maxunit::Int
  verbose::Bool
  verbosestep::Int
end

function KSettings(outfile::AbstractString,
                   T::Int,
                   burnin::Int,
                   tstep::Int,
                   op::Vector{Float64},
                   α::Float64,
                   θ::Float64,
                   γ::Vector{Float64},
                   r::Float64,
                   λs1::Float64,
                   λs2::Float64,
                   parawm::Float64,
                   maxclust::Int,
                   maxunit::Int,
                   verbose::Bool,
                   verbosestep::Int)
  KSettings(outfile, T, burnin, tstep, StatsBase.WeightVec(op), α, θ, γ, r,
            Distributions.Beta(λs1, λs2), parawm, maxclust, maxunit, verbose,
            verbosestep)
end
