# This file is part of Kpax3. License is MIT.

function kpax3(x::AminoAcidData,
               partition::AbstractString,
               settings::KSettings)
  R = normalizepartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(k, settings.maxclust))

  support = KSupport(size(x.data, 1), size(x.data, 2), settings.maxclust,
                     settings.maxunit)

  state = AminoAcidState(x.data, R, priorR, priorC, settings)

  kpax3mcmc!(x.data, priorR, priorC, settings, support, state)

  nothing
end
