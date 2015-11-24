# This file is part of Kpax3. License is MIT.

function kpax3aa(x::AminoAcidData,
                 partition::AbstractString,
                 T::Int,
                 outfile::AbstractString;
                 burnin::Int=0,
                 op::Array{Float64, 1}=[1.0, 0.0, 0.0, 0.0],
                 t::Int=1,
                 maxclust::Int=500,
                 α::Float64=0.0,
                 θ::Float64=1.0,
                 γ::Array{Float64, 1}=[0.6, 0.35, 0.05],
                 r::Float64=log(0.001) / log(0.95),
                 verbose::Bool=true,
                 verbosestep::Int=100000)
  R = initpartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(α, θ)
  priorC = AminoAcidPriorCol(x.data, k, γ, r)

  # To preallocate memory, we must choose an upper bound to the maximum
  # number of clusters. The total number of clusters can't be greater than n,
  # but n could be a very very large number. Choose maxclust as the upper
  # bound unless the current partition has more than maxclust clusters...
  maxclust = max(min(maxclust, n), k)

  aamcmc = AminoAcidMCMC(x.data, R, priorC, T, maxclust, burnin, op, t, outfile)

  kpax3mcmc(aamcmc, priorR, priorC, x.data, verbose, verbosestep)
end
