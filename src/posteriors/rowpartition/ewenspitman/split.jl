# This file is part of Kpax3. License is MIT.

function logSplitRatioPriorRow(ep::EwensPitmanPAUT,
                               v::Array{Float64, 1},
                               k::Float64)
  log(ep.θ + k * ep.α) + lgamma(v[1] - ep.α) + lgamma(v[2] - ep.α) -
  log(1.0 - ep.α) - lgamma(v[3] - ep.α)
end

function logSplitRatioPriorRow(ep::EwensPitmanPAZT,
                               v::Array{Float64, 1},
                               k::Float64)
  log(k * ep.α) + lgamma(v[1] - ep.α) + lgamma(v[2] - ep.α) - log(1.0 - ep.α) -
  lgamma(v[3] - ep.α)
end

function logSplitRatioPriorRow(ep::EwensPitmanZAPT,
                               v::Array{Float64, 1},
                               k::Float64)
  log(ep.θ) + lgamma(v[1]) + lgamma(v[2]) - lgamma(v[3])
end

# TODO: risk of overflow
function logSplitRatioPriorRow(ep::EwensPitmanPAUT,
                               v::Array{Float64, 1},
                               k::Float64)
  if v[1] > v[2]
    log((k - ep.L) * ep.α * prod(1:(v[2] - 1) - ep.α) /
        prod(v[1]:(v[3] - 1) - ep.α))
  else
    log((k - ep.L) * ep.α * prod(1:(v[1] - 1) - ep.α) /
        prod(v[2]:(v[3] - 1) - ep.α))
  end
end
