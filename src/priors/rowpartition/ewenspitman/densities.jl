# This file is part of Kpax3. License is MIT.

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
function dPriorRow(ep::EwensPitmanPAUT,
                   p::Array{Int, 1})
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

function logdPriorRow(ep::EwensPitmanPAUT,
                      p::Array{Int, 1})
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

function dPriorRow(ep::EwensPitmanPAUT,
                   n::Int,
                   k::Int,
                   m::Array{Int, 1})
  exp((k - 1) * log(ep.α) + lgamma(ep.θ / ep.α + k) - lgamma(ep.θ / ep.α + 1) +
      lgamma(ep.θ + 1) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0] - ep.α)) -
      k * lgamma(1 - ep.α))
end

function logdPriorRow(ep::EwensPitmanPAUT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  (k - 1) * log(ep.α) + lgamma(ep.θ / ep.α + k) - lgamma(ep.θ / ep.α + 1) +
  lgamma(ep.θ + 1) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0] - ep.α)) -
  k * lgamma(1 - ep.α)
end

function dPriorRow(ep::EwensPitmanPAZT,
                   p::Array{Int, 1})
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

function logdPriorRow(ep::EwensPitmanPAZT,
                      p::Array{Int, 1})
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

function dPriorRow(ep::EwensPitmanPAZT,
                   n::Int,
                   k::Int,
                   m::Array{Int, 1})
  exp((k - 1) * log(ep.α) + lgamma(k) - lgamma(n) +
      sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α))
end

function logdPriorRow(ep::EwensPitmanPAZT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  (k - 1) * log(ep.α) + lgamma(k) - lgamma(n) + sum(lgamma(m[m .> 0] - ep.α)) -
  k * lgamma(1 - ep.α)
end

function dPriorRow(ep::EwensPitmanZAPT,
                   p::Array{Int, 1})
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

function logdPriorRow(ep::EwensPitmanZAPT,
                      p::Array{Int, 1})
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

function dPriorRow(ep::EwensPitmanZAPT,
                   n::Int,
                   k::Int,
                   m::Array{Int, 1})
  exp(k * log(ep.θ) + lgamma(ep.θ) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0])))
end

function logdPriorRow(ep::EwensPitmanZAPT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  k * log(ep.θ) + lgamma(ep.θ) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0]))
end

# TODO: risk of integer overflow
# Is it possible to avoid it without introducing rounding errors or numerical
# instability?
function dPriorRow(ep::EwensPitmanNAPT,
                   p::Array{Int, 1})
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

function logdPriorRow(ep::EwensPitmanNAPT,
                      p::Array{Int, 1})
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

function dPriorRow(ep::EwensPitmanNAPT,
                   n::Int,
                   k::Int,
                   m::Array{Int, 1})
  prod((1:(k - 1)) - ep.L) * ep.α^(k - 1) *
  exp(sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α)) /
  prod((1:(n - 1)) - ep.α * ep.L)
end

function logdPriorRow(ep::EwensPitmanNAPT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  log(prod((1:(k - 1)) - ep.L) * ep.α^(k - 1) *
      exp(sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α)) /
      prod((1:(n - 1)) - ep.α * ep.L))
end
