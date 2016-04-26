# This file is part of Kpax3. License is MIT.

"""
# User defined settings for a Kpax3 run

## Description

## Fields

"""
immutable KSettings
  fpath::AbstractString
  T::Int
  burnin::Int
  tstep::Int
  op::StatsBase.WeightVec
  α::Real
  θ::Real
  γ::Vector{Float64}
  r::Float64
  distws::Distributions.Beta
  parawm::Float64
  maxclust::Int
  maxunit::Int
  verbose::Bool
  verbosestep::Int
end

function KSettings(fpath::AbstractString;
                   T::Int=10000,
                   burnin::Int=1000,
                   tstep::Int=1,
                   op::Vector{Float64}=[0.6; 0.3; 0.1],
                   α::Real=0.5,
                   θ::Real=-0.1,
                   γ::Vector{Float64}=[0.6; 0.35; 0.05],
                   r::Float64=log(0.001) / log(0.95),
                   λs1::Float64=1.0,
                   λs2::Float64=1.0,
                   parawm::Float64=5.0,
                   maxclust::Int=500,
                   maxunit::Int=500,
                   verbose::Bool=true,
                   verbosestep::Int=500)
  # open oufile for writing and immediately close it. We do this to throw a
  # proper Julia standard exception if something is wrong
  f = open(fpath, "a")
  close(f)

  if T < 1
    throw(KDomainError("Argument 'T' is lesser than 1."))
  end

  if burnin < 0
    throw(KDomainError("Argument 'burnin' is negative."))
  end

  if tstep < 0
    throw(KDomainError("Argument 'tstep' is negative."))
  end

  if length(op) != 3
    throw(KInputError("Argument 'op' does not have length 3."))
  elseif op[1] < 0
    throw(KDomainError("Argument 'op[1]' is negative."))
  elseif op[2] < 0
    throw(KDomainError("Argument 'op[2]' is negative."))
  elseif op[3] < 0
    throw(KDomainError("Argument 'op[3]' is negative."))
  end

  if length(γ) != 3
    throw(KInputError("Argument 'γ' does not have length 3."))
  elseif γ[1] < 0
    throw(KDomainError("Argument 'γ[1]' is negative."))
  elseif γ[2] < 0
    throw(KDomainError("Argument 'γ[2]' is negative."))
  elseif γ[3] < 0
    throw(KDomainError("Argument 'γ[3]' is negative."))
  end

  if r <= 0.0
    throw(KDomainError("Argument 'r' is not positive."))
  end

  if λs1 <= 0.0
    throw(KDomainError("Argument 'λs1' is not positive."))
  end

  if λs2 <= 0.0
    throw(KDomainError("Argument 'λs2' is not positive."))
  end

  if parawm <= 0.0
    throw(KDomainError("Argument 'parawm' is not positive."))
  end

  if maxclust < 1
    throw(KDomainError("Argument maxclust is lesser than 1."))
  end

  if maxunit < 1
    throw(KDomainError("Argument maxunit is lesser than 1."))
  end

  if verbosestep < 0
    throw(KDomainError("Argument 'verbosestep' is negative."))
  end

  KSettings(fpath, T, burnin, tstep, StatsBase.WeightVec(op), α, θ, γ, r,
            Distributions.Beta(λs1, λs2), parawm, maxclust, maxunit, verbose,
            verbosestep)
end
