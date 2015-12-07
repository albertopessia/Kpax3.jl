# This file is part of Kpax3. License is MIT.

type KWeight
  w::Array{Float64, 2}
  c::Array{Float64, 1}
  z::Array{Float64, 2}
end

"""
# Support object

## Description

An object of type KSupport contains data structures (matrices, vectors, etc.)
used during numerical computations. We are trying to avoid dynamic memory
allocation of the same objects over and over again.

"""
type KSupport
  gi::KCluster
  wi::KWeight

  gj::KCluster
  wj::KWeight
end

function KSupport(m::Int,
                  n::Int)
  gi = KCluster(zero(Int), zeros(Int, 1), zeros(Float64, m))
  wi = KWeight(zeros(Float64, 4, m), zeros(Float64, m), zeros(Float64, 4, m))

  gj = KCluster(zero(Int), zeros(Int, 1), zeros(Float64, m))
  wj = KWeight(zeros(Float64, 4, m), zeros(Float64, m), zeros(Float64, 4, m))

  KSupport(gi, wi, gj, wj)
end
