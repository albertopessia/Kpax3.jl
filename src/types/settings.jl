# This file is part of Kpax3. License is MIT.

"""
# User defined settings for a Kpax3 run

## Description

## Fields

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
