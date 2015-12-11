# This file is part of Kpax3. License is MIT.

function kpax3aa(x::AminoAcidData,
                 partition::AbstractString,
                 outfile::AbstractString,
                 T::Int;
                 burnin::Int=0,
                 tstep::Int=1,
                 op::Vector{Float64}=[1.0; 0.0; 0.0],
                 α::Float64=0.0,
                 θ::Float64=1.0,
                 γ::Vector{Float64}=[0.6; 0.35; 0.05],
                 r::Float64=log(0.001) / log(0.95),
                 λs1::Float64=1.0,
                 λs2::Float64=1.0,
                 parawm::Float64=5.0,
                 maxclust::Int=500,
                 maxunit::Int=500,
                 verbose::Bool=true,
                 verbosestep::Int=100000)
  R = initpartition(partition, x.id)
  k = maximum(R)

  # To preallocate memory, we must choose an upper bound to the maximum
  # number of clusters. The total number of clusters can't be greater than n,
  # but n could be a very very large number. Choose maxclust as the upper
  # bound unless the current partition has more than maxclust clusters...
  maxclust = max(min(maxclust, size(x.data, 2)), k)
  maxunit = min(maxunit, size(x.data, 2))

  settings = KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r, λs1, λs2,
                       parawm, maxclust, maxunit, verbose, verbosestep)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, k, settings.γ, settings.r)

  support = KSupport(size(x.data, 1), size(x.data, 2), settings.maxclust,
                     settings.maxunit)

  mcmcobj = AminoAcidMCMC(x.data, R, priorR, priorC, settings)

  kpax3mcmc!(x.data, priorR, priorC, settings, support, mcmcobj)
end
