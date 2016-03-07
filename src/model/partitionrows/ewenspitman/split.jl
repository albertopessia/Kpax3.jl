# This file is part of Kpax3. License is MIT.

function logratiopriorrowsplit(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAUT)
  log(ep.θ + (k - 1) * ep.α) - lgamma(1 - ep.α) +
  lgamma(vi - ep.α) + lgamma(vj - ep.α) - lgamma(vi + vj - ep.α)
end

function logratiopriorrowsplit(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAZT)
  log(k - 1) + log(ep.α) - lgamma(1 - ep.α) +
  lgamma(vi - ep.α) + lgamma(vj - ep.α) - lgamma(vi + vj - ep.α)
end

function logratiopriorrowsplit(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanZAPT)
  log(ep.θ) + lgamma(vi) + lgamma(vj) - lgamma(vi + vj)
end

function logratiopriorrowsplit(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanNAPT)
  log((k - 1 - ep.L) * ep.α * exp(lgamma(vi - ep.α) + lgamma(vj - ep.α) -
      lgamma(1 - ep.α) - lgamma(vi + vj - ep.α)))
end
