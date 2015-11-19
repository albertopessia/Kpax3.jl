"""
# Density of the Ewens-Pitman distribution

## Description

Probability of a partition according to the Ewens-Pitman distribution.

## Usage

dEwensPitman(ep, p)
dEwensPitman(ep, n, k, m)

## Arguments

* `ep` Object of (super)type EwensPitman
* `p` Vector of integers representing a partition
* `n` Set size (Integer)
* `k` Number of blocks (Integer)
* `m` Vector of integers representing block sizes

## Details

## Examples

"""
function dEwensPitman(ep::EwensPitmanPAUT,
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

function dEwensPitman(ep::EwensPitmanPAUT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  exp((k - 1) * log(ep.α) + lgamma(ep.θ / ep.α + k) - lgamma(ep.θ / ep.α + 1) +
      lgamma(ep.θ + 1) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0] - ep.α)) -
      k * lgamma(1 - ep.α))
end

function dEwensPitman(ep::EwensPitmanPAZT,
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

function dEwensPitman(ep::EwensPitmanPAZT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  exp((k - 1) * log(ep.α) + lgamma(k) - lgamma(n) +
      sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α))
end

function dEwensPitman(ep::EwensPitmanZAPT,
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

function dEwensPitman(ep::EwensPitmanZAPT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  exp(k * log(ep.θ) + lgamma(ep.θ) - lgamma(ep.θ + n) + sum(lgamma(m[m .> 0])))
end

# TODO: risk of integer overflow
# Is it possible to avoid it without introducing rounding errors or numerical
# instability?
function dEwensPitman(ep::EwensPitmanNAPT,
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

function dEwensPitman(ep::EwensPitmanNAPT,
                      n::Int,
                      k::Int,
                      m::Array{Int, 1})
  prod((1:(k - 1)) - ep.L) * ep.α^(k - 1) *
  exp(sum(lgamma(m[m .> 0] - ep.α)) - k * lgamma(1 - ep.α)) /
  prod((1:(n - 1)) - ep.α * ep.L)
end
