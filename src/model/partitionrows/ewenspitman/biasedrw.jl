# This file is part of Kpax3. License is MIT.

function logratiopriorrowbrwmerge(k::Real,
                                  vj::Real,
                                  ep::EwensPitmanPAUT)
  log(vj - ep.α) - log(ep.θ + k * ep.α)
end

function logratiopriorrowbrwmove(vi::Real,
                                 vj::Real,
                                 ep::EwensPitmanPAUT)
  log(vj - ep.α) - log(vi - 1 - ep.α)
end

function logratiopriorrowbrwsplit(k::Real,
                                  vi::Real,
                                  ep::EwensPitmanPAUT)
  log(ep.θ + (k - 1) * ep.α) - log(vi - 1 - ep.α)
end

function logratiopriorrowbrwmerge(k::Real,
                                  vj::Real,
                                  ep::EwensPitmanPAZT)
  log(vj - ep.α) - log(k) - log(ep.α)
end

function logratiopriorrowbrwmove(vi::Real,
                                 vj::Real,
                                 ep::EwensPitmanPAZT)
  log(vj - ep.α) - log(vi - 1 - ep.α)
end

function logratiopriorrowbrwsplit(k::Real,
                                  vi::Real,
                                  ep::EwensPitmanPAZT)
  log(k - 1) + log(ep.α) - log(vi - 1 - ep.α)
end

function logratiopriorrowbrwmerge(k::Real,
                                  vj::Real,
                                  ep::EwensPitmanZAPT)
  log(vj) - log(ep.θ)
end

function logratiopriorrowbrwmove(vi::Real,
                                 vj::Real,
                                 ep::EwensPitmanZAPT)
  log(vj) - log(vi - 1)
end

function logratiopriorrowbrwsplit(k::Real,
                                  vi::Real,
                                  ep::EwensPitmanZAPT)
  log(ep.θ) - log(vi - 1)
end

function logratiopriorrowbrwmerge(k::Real,
                                  vj::Real,
                                  ep::EwensPitmanNAPT)
  - log((k - ep.L) * ep.α * exp(lgamma(vj - ep.α) - lgamma(vj + 1 - ep.α)))
end

function logratiopriorrowbrwmove(vi::Real,
                                 vj::Real,
                                 ep::EwensPitmanNAPT)
  log(vj - ep.α) - log(vi - 1 - ep.α)
end

function logratiopriorrowbrwsplit(k::Real,
                                  vi::Real,
                                  ep::EwensPitmanNAPT)
  log((k - 1 - ep.L) * ep.α * exp(lgamma(vi - 1 - ep.α) - lgamma(vi - ep.α)))
end
