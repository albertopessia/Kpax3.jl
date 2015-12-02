# This file is part of Kpax3. License is MIT.

function logdPriorRow(p::Array{Int, 1},
                      ep::EwensPitmanPAUT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  (k - 1) * log(ep.α) + lgamma(ep.θ / ep.α + k) - lgamma(ep.θ / ep.α + 1) +
  lgamma(ep.θ + 1) - lgamma(ep.θ + n) + sum(lgamma(m[idx] - ep.α)) -
  k * lgamma(1 - ep.α)
end

function logdPriorRow(n::Int,
                      k::Int,
                      m::Array{Int, 1},
                      ep::EwensPitmanPAUT)
  (k - 1) * log(ep.α) + lgamma(ep.θ / ep.α + k) - lgamma(ep.θ / ep.α + 1) +
  lgamma(ep.θ + 1) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0] - ep.α)) -
  k * lgamma(1 - ep.α)
end

function logdPriorRow(p::Array{Int, 1},
                      ep::EwensPitmanPAZT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  (k - 1) * log(ep.α) + lgamma(k) - lgamma(n) + sum(lgamma(m[idx] - ep.α)) -
  k * lgamma(1 - ep.α)
end

function logdPriorRow(n::Int,
                      k::Int,
                      m::Array{Int, 1},
                      ep::EwensPitmanPAZT)
  (k - 1) * log(ep.α) + lgamma(k) - lgamma(n) + sum(lgamma(m[m .> 0] - ep.α)) -
  k * lgamma(1 - ep.α)
end

function logdPriorRow(p::Array{Int, 1},
                      ep::EwensPitmanZAPT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  k * log(ep.θ) + lgamma(ep.θ) - lgamma(ep.θ + n) + sum(lgamma(m[idx]))
end

function logdPriorRow(n::Int,
                      k::Int,
                      m::Array{Int, 1},
                      ep::EwensPitmanZAPT)
  k * log(ep.θ) + lgamma(ep.θ) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0]))
end

function logdPriorRow(p::Array{Int, 1},
                      ep::EwensPitmanNAPT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  log(prod((1:(k - 1)) - ep.L) * ep.α^(k - 1) *
      exp(sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α)) /
      prod((1:(n - 1)) - ep.α * ep.L))
end

function logdPriorRow(n::Int,
                      k::Int,
                      m::Array{Int, 1},
                      ep::EwensPitmanNAPT)
  log(prod((1:(k - 1)) - ep.L) * ep.α^(k - 1) *
      exp(sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α)) /
      prod((1:(n - 1)) - ep.α * ep.L))
end

"""
# Density of the Ewens-Pitman distribution

## Description

Probability of a partition according to the Ewens-Pitman distribution.

## Usage

dPriorRow(ep, p)
dPriorRow(ep, n, k, m)

## Arguments

* `ep` Object of (super)type EwensPitman
* `p` Vector of integers representing a partition
* `n` Set size (Integer)
* `k` Number of blocks (Integer)
* `m` Vector of integers representing block sizes

## Details

## Examples

"""
function dPriorRow(p::Array{Int, 1},
                   ep::EwensPitmanPAUT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  exp((k - 1) * log(ep.α) + lgamma(ep.θ / ep.α + k) - lgamma(ep.θ / ep.α + 1) +
      lgamma(ep.θ + 1) - lgamma(ep.θ + n) + sum(lgamma(m[idx] - ep.α)) -
      k * lgamma(1 - ep.α))
end

function dPriorRow(n::Int,
                   k::Int,
                   m::Array{Int, 1},
                   ep::EwensPitmanPAUT)
  exp((k - 1) * log(ep.α) + lgamma(ep.θ / ep.α + k) - lgamma(ep.θ / ep.α + 1) +
      lgamma(ep.θ + 1) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0] - ep.α)) -
      k * lgamma(1 - ep.α))
end

function dPriorRow(p::Array{Int, 1},
                   ep::EwensPitmanPAZT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  exp((k - 1) * log(ep.α) + lgamma(k) - lgamma(n) + sum(lgamma(m[idx] - ep.α)) -
      k * lgamma(1 - ep.α))
end

function dPriorRow(n::Int,
                   k::Int,
                   m::Array{Int, 1},
                   ep::EwensPitmanPAZT)
  exp((k - 1) * log(ep.α) + lgamma(k) - lgamma(n) +
      sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α))
end

function dPriorRow(p::Array{Int, 1},
                   ep::EwensPitmanZAPT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  exp(k * log(ep.θ) + lgamma(ep.θ) - lgamma(ep.θ + n) + sum(lgamma(m[idx])))
end

function dPriorRow(n::Int,
                   k::Int,
                   m::Array{Int, 1},
                   ep::EwensPitmanZAPT)
  exp(k * log(ep.θ) + lgamma(ep.θ) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0])))
end

# TODO: risk of integer overflow
# Is it possible to avoid it without introducing rounding errors or numerical
# instability?
function dPriorRow(p::Array{Int, 1},
                   ep::EwensPitmanNAPT)
  n = length(p)
  k = 0.0

  m = zeros(Float64, n)
  idx = falses(n)

  for i in 1:n
    m[p[i]] += 1.0

    if m[p[i]] == 1.0
      idx[p[i]] = true
      k += 1.0
    end
  end

  prod((1:(k - 1)) - ep.L) * ep.α^(k - 1) *
  exp(sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α)) /
  prod((1:(n - 1)) - ep.α * ep.L)
end

function dPriorRow(n::Int,
                   k::Int,
                   m::Array{Int, 1},
                   ep::EwensPitmanNAPT)
  prod((1:(k - 1)) - ep.L) * ep.α^(k - 1) *
  exp(sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α)) /
  prod((1:(n - 1)) - ep.α * ep.L)
end
