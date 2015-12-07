# This file is part of Kpax3. License is MIT.

function logMergeRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAUT)
  log(1.0 - ep.α) + lgamma(vi + vj - ep.α) - log(ep.θ + k * ep.α) -
  lgamma(vi - ep.α) - lgamma(vj - ep.α)
end

function logMergeRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAZT)
  log(1.0 - ep.α) + lgamma(vi + vj - ep.α) - log(k * ep.α) - lgamma(vi - ep.α) -
  lgamma(vj - ep.α)
end

function logMergeRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanZAPT)
  lgamma(vi + vj) - log(ep.θ) - lgamma(vi) - lgamma(vj)
end

# TODO: risk of overflow
function logMergeRatioPriorRow(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAUT)
  v = vi + vj

  if vi > vj
    log(prod(vi:(v - 1) - ep.α) / ((k - ep.L) * ep.α * prod(1:(vj - 1) - ep.α)))
  else
    log(prod(vj:(v - 1) - ep.α) / ((k - ep.L) * ep.α * prod(1:(vi - 1) - ep.α)))
  end
end
