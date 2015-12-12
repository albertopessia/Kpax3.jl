# This file is part of Kpax3. License is MIT.

type KWeight
  c::Vector{Float64}
  w::Matrix{Float64}
  z::Matrix{Float64}
end

"""
# Support object

## Description

An object of type KSupport contains data structures (matrices, vectors, etc.)
used during numerical computations. We are trying to avoid dynamic memory
allocation of the same objects over and over again.

"""
type KSupport
  vi::Int
  ni::Vector{Float64}
  ui::Vector{Int}
  wi::KWeight

  vj::Int
  nj::Vector{Float64}
  uj::Vector{Int}
  wj::KWeight

  C::Matrix{UInt8}
end

function KSupport(m::Int,
                  n::Int,
                  maxclust::Int,
                  maxunit::Int)
  vi = 0
  ni = zeros(Float64, m)
  ui = zeros(Int, maxunit)
  wi = KWeight(zeros(Float64, m), zeros(Float64, 4, m), zeros(Float64, 4, m))

  vj = 0
  nj = zeros(Float64, m)
  uj = zeros(Int, maxunit)
  wj = KWeight(zeros(Float64, m), zeros(Float64, 4, m), zeros(Float64, 4, m))

  C = zeros(UInt8, maxclust, m)

  KSupport(vi, ni, ui, wi, vj, nj, uj, wj, C)
end
