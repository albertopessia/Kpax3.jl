# This file is part of Kpax3. License is MIT.

function kpax3aa(x::AminoAcidData,
                 partition::AbstractString,
                 T::Int,
                 outfile::AbstractString;
                 burnin::Int=0,
                 t::Int=1,
                 op::Array{Float64, 1}=[0.5, 0.25, 0.25],
                 α::Float64=0.0,
                 θ::Float64=1.0,
                 γ::Array{Float64, 1}=[0.6, 0.35, 0.05],
                 r::Float64=log(0.001) / log(0.95),
                 λs1::Float64=1.0,
                 λs2::Float64=1.0,
                 parawm::Float64=5.0,
                 maxclust::Int=500,
                 maxunit::Int=500,
                 verbose::Bool=true,
                 verbosestep::Int=100000)
  settings = Kpax3Settings(T, outfile, burnin, t, op, α, θ, γ, r, λs1, λs2,
                           parawm, maxclust, maxunit, verbose, verbosestep)

  R = initpartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, k, settings.γ, settings.r)

  # To preallocate memory, we must choose an upper bound to the maximum
  # number of clusters. The total number of clusters can't be greater than n,
  # but n could be a very very large number. Choose maxclust as the upper
  # bound unless the current partition has more than maxclust clusters...
  maxclust = max(min(maxclust, size(x.data, 2)), k)

  mcmcobj = AminoAcidMCMC(x.data, R, priorC, settings)

  kpax3mcmc(mcmcobj, priorR, priorC, x.data, settings)
end
