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
  m::Int
  n::Int

  vi::Int
  ni::Vector{Float64}
  ui::Vector{Int}
  wi::KWeight

  vj::Int
  nj::Vector{Float64}
  uj::Vector{Int}
  wj::KWeight

  tmp::Vector{Float64}

  cl::Vector{Int}
  k::Int

  C::Matrix{UInt8}

  logω::Vector{Float64}

  logpC::Vector{Float64}
  lograR::Float64
  loglik::Float64
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

  tmp = zeros(Float64, 4)

  cl = zeros(Int, n)
  k = 0

  C = zeros(UInt8, maxclust, m)

  logω = zeros(Float64, 4)

  logpC = zeros(Float64, 2)

  KSupport(m, n, vi, ni, ui, wi, vj, nj, uj, wj, tmp, cl, k, C, logω, logpC,
           0.0, 0.0)
end
