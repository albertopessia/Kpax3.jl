# This file is part of Kpax3. License is MIT.

"""
# User defined settings for a Kpax3 run

## Description

## Fields

"""
immutable KSettings
  ifile::AbstractString
  ofile::AbstractString
  # Common parameters
  protein::Bool
  miss::Vector{UInt8}
  l::Int
  α::Real
  θ::Real
  γ::Vector{Float64}
  r::Float64
  maxclust::Int
  maxunit::Int
  verbose::Bool
  verbosestep::Int
  # Genetic Algorithm parameters
  popsize::Int
  maxiter::Int
  maxgap::Int
  xrate::Float64
  mrate::Float64
  # MCMC parameters
  T::Int
  burnin::Int
  tstep::Int
  op::StatsBase.WeightVec
  distws::Distributions.Beta
  parawm::Float64
end

function KSettings(ifile::AbstractString,
                   ofile::AbstractString;
                   protein::Bool=true,
                   miss::Vector{UInt8}=zeros(UInt8, 0),
                   l::Int=100000000,
                   α::Real=0.5,
                   θ::Real=-0.1,
                   γ::Vector{Float64}=[0.6; 0.35; 0.05],
                   r::Float64=log(0.001) / log(0.95),
                   maxclust::Int=500,
                   maxunit::Int=500,
                   verbose::Bool=false,
                   verbosestep::Int=500,
                   popsize::Int=50,
                   maxiter::Int=10000,
                   maxgap::Int=10000,
                   xrate::Float64=0.9,
                   mrate::Float64=0.005,
                   T::Int=1000000,
                   burnin::Int=100000,
                   tstep::Int=1,
                   op::Vector{Float64}=[0.6; 0.3; 0.1],
                   λs1::Float64=1.0,
                   λs2::Float64=1.0,
                   parawm::Float64=5.0)
  # open files and immediately close them. We do this to throw a proper Julia
  # standard exception if something is wrong
  f = open(ifile, "r")
  close(f)

  f = open(ofile, "a")
  close(f)

  if length(miss) == 0
    miss = if protein
             UInt8['?', '*', '#', '-', 'b', 'j', 'x', 'z']
           else
             UInt8['?', '*', '#', '-', 'b', 'd', 'h', 'k', 'm', 'n', 'r', 's',
                   'v', 'w', 'x', 'y', 'j', 'z']
           end
  else
    if length(miss) == 1
      if miss[1] != UInt8(0)
        if UInt8(0) < miss[1] < UInt8(128)
          if UInt8(64) < miss[1] < UInt8(91)
            # convert to lowercase
            miss[1] += UInt8(32)
          end
        else
          throw(KDomainError(string("Value 'miss[1]' is not in the range ",
                                    "[1, ..., 127]: ", Int(miss[1]), ".")))
        end
      end
    else
      for i in 1:length(miss)
        if UInt8(0) < miss[i] < UInt8(128)
          if UInt8(64) < miss[i] < UInt8(91)
            # convert to lowercase
            miss[i] += UInt8(32)
          end
        else
          throw(KDomainError(string("Value 'miss[", i, "]' is not in the ",
                                    "range [1, ..., 127]: ",
                                    Int(miss[i]), ".")))
        end
      end
    end
  end

  if l < 1
    throw(KDomainError("Argument 'l' is not positive."))
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

  if maxclust < 1
    throw(KDomainError("Argument maxclust is lesser than 1."))
  end

  if maxunit < 1
    throw(KDomainError("Argument maxunit is lesser than 1."))
  end

  if verbosestep < 0
    throw(KDomainError("Argument 'verbosestep' is negative."))
  end

  # disable status reports if verbosestep is not positive
  verbose = verbose && (verbosestep > 0)

  if popsize < 4
    throw(KDomainError("Argument 'popsize' is lesser than 4."))
  end

  if maxiter < 1
    throw(KDomainError("Argument 'maxiter' is lesser than 1."))
  end

  if maxgap < 0
    throw(KDomainError("Argument 'maxgap' is lesser than 1."))
  end

  if !(0 <= xrate <= 1)
    throw(KDomainError("Argument 'xrate' is not in the range [0, 1]."))
  end

  if !(0 <= mrate <= 1)
    throw(KDomainError("Argument 'mrate' is not in the range [0, 1]."))
  end

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

  if λs1 <= 0.0
    throw(KDomainError("Argument 'λs1' is not positive."))
  end

  if λs2 <= 0.0
    throw(KDomainError("Argument 'λs2' is not positive."))
  end

  if parawm <= 0.0
    throw(KDomainError("Argument 'parawm' is not positive."))
  end

  KSettings(ifile, ofile, protein, miss, l, α, θ, γ, r, maxclust, maxunit,
            verbose, verbosestep, popsize, maxiter, maxgap, xrate, mrate, T,
            burnin, tstep, StatsBase.WeightVec(op),
            Distributions.Beta(λs1, λs2), parawm)
end
