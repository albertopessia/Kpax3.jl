# This file is part of Kpax3. License is MIT.

"""
# Support object

## Description

An object of type Kpax3Support contains data structures (matrices, vectors,
etc.) used during numerical computations. We are trying to avoid dynamic
memory allocation of the same objects over and over again.

"""
type Kpax3Support
  smv::Array{Float64, 1}
  smR::Array{Int, 1}
  smMPi::Array{Float64, 2}
  smMPj::Array{Float64, 2}
  smCP::Array{Float64, 1}
end

function Kpax3Support(m::Int,
                      n::Int)
  smv = zeros(Float64, 3)
  smR = zeros(Int, n)

  smMPi = zeros(Float64, m, 4)
  smMPj = zeros(Float64, m, 4)

  smCP = zeros(Float64, n)

  Kpax3Support(smv, smR, smMPi, smMPj, smCP)
end
