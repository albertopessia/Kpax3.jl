# This file is part of Kpax3. License is MIT.

function logratiopriorrowmerge(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAUT)
  lgamma(1 - ep.α) + lgamma(vi + vj - ep.α) - log(ep.θ + k * ep.α) -
  lgamma(vi - ep.α) - lgamma(vj - ep.α)
end

function logratiopriorrowmerge(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanPAZT)
  lgamma(1 - ep.α) + lgamma(vi + vj - ep.α) - log(k * ep.α) - lgamma(vi - ep.α) -
  lgamma(vj - ep.α)
end

function logratiopriorrowmerge(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanZAPT)
  lgamma(vi + vj) - log(ep.θ) - lgamma(vi) - lgamma(vj)
end

function logratiopriorrowmerge(k::Real,
                               vi::Real,
                               vj::Real,
                               ep::EwensPitmanNAPT)
  - log((k - ep.L) * ep.α * exp(lgamma(vi - ep.α) + lgamma(vj - ep.α) -
        lgamma(1 - ep.α) - lgamma(vi + vj - ep.α)))
end
