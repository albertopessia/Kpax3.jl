# This file is part of Kpax3. License is MIT.

function logSplitRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAUT)
  log(ep.θ + k * ep.α) + lgamma(vi - ep.α) + lgamma(vj - ep.α) -
  log(1.0 - ep.α) - lgamma(vi + vj - ep.α)
end

function logSplitRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAZT)
  log(k * ep.α) + lgamma(vi - ep.α) + lgamma(vj - ep.α) - log(1.0 - ep.α) -
  lgamma(vi + vj - ep.α)
end

function logSplitRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanZAPT)
  log(ep.θ) + lgamma(vi) + lgamma(vj) - lgamma(vi + vj)
end

# TODO: risk of overflow
function logSplitRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAUT)
  v = vi + vj

  if vi > vj
    log((k - ep.L) * ep.α * prod(1:(vj - 1) - ep.α) / prod(vi:(v - 1) - ep.α))
  else
    log((k - ep.L) * ep.α * prod(1:(vi - 1) - ep.α) / prod(vj:(v - 1) - ep.α))
  end
end
