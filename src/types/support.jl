# This file is part of Kpax3. License is MIT.

"""
# Support object

## Description

An object of type Kpax3Support contains data structures (matrices, vectors,
etc.) used during numerical computations. We are trying to avoid dynamic
memory allocation of the same objects over and over again.

"""
type Kpax3Support
  smR::Array{Float64, 1}
  smV::Array{Float64, 1}
  smCO::Array{Float64, 2}
  smMP::Array{Float64, 3}
  smLC::Array{Float64, 2}
  smCP::Array{Float64, 1}
end

function Kpax3Support(m::Int,
                      n::Int)
  smR = zeros(Float64, n)
  smV = zeros(Float64, 3)
  smCO = zeros(Float64, 2, m)
  smMP = zeros(Float64, 4, 2, m)
  smLC = zeros(Float64, 2, m)
  smCP = zeros(Float64, n)

  Kpax3Support(smR, smV, smCO, smMP, smLC, smCP)
end
