# This file is part of Kpax3. License is MIT.

function split!(mcmcobj::AminoAcidMCMC,
                priorR::PriorRowPartition,
                priorC::PriorColPartition,
                data::Array{UInt8, 2},
                ij::Array{Int, 1},
                neighbours::Array{Int, 1},
                S::Int,
                settings::Kpax3Settings,
                support::Kpax3Support)
  # initialize variables
  gi = mcmcobj.R[ij[1]]

  # sample sizes and partition allocations
  support.smv[1] = one(Float64)
  support.smv[2] = one(Float64)
  support.smv[3] = S + 2.0

  fill!(support.smR, zero(Int))

  fill!(support.smMPi, zero(Float64))
  fill!(support.smMPj, zero(Float64))

  fill!(support.smCP, zero(Float64))

  for j in 1:m
    support.smMPi[j, 1] = exp(priorC.logγ[1] +
                              logmarglik(data[j, ij[1]], 1, priorC.A[j, 1],
                                         priorC.B[j, 1]))

    support.smMPi[j, 2] = exp(priorC.logγ[2] +
                              logmarglik(data[j, ij[1]], 1, priorC.A[j, 2],
                                         priorC.B[j, 2]))

    support.smMPi[j, 3] = exp(priorC.logγ[3] +
                              log(mcmcobj.k) - log(mcmcobj.k + 1) +
                              logmarglik(data[j, ij[1]], 1, priorC.A[j, 3],
                                         priorC.B[j, 3]))

    support.smMPi[j, 4] = exp(priorC.logγ[3] -
                              log(mcmcobj.k + 1) +
                              logmarglik(data[j, ij[1]], 1, priorC.A[j, 4],
                                         priorC.B[j, 4]))

    support.smMPj[j, 1] = exp(priorC.logγ[1] +
                              logmarglik(data[j, ij[2]], 1, priorC.A[j, 1],
                                         priorC.B[j, 1]))

    support.smMPj[j, 2] = exp(priorC.logγ[2] +
                              logmarglik(data[j, ij[2]], 1, priorC.A[j, 2],
                                         priorC.B[j, 2]))

    support.smMPj[j, 3] = exp(priorC.logγ[3] +
                              log(mcmcobj.k) - log(mcmcobj.k + 1) +
                              logmarglik(data[j, ij[2]], 1, priorC.A[j, 3],
                                         priorC.B[j, 3]))

    support.smMPj[j, 4] = exp(priorC.logγ[3] -
                              log(mcmcobj.k + 1) +
                              logmarglik(data[j, ij[2]], 1, priorC.A[j, 4],
                                         priorC.B[j, 4]))
  end

  # log-likelihoods
  ll = zeros(Float64, 3)

  # logarithm of the sequential allocation probabilities
  lq = zero(Float64)

  # support weight vectors
  wv1 = zeros(Float64, 2)
  wv2 = zeros(Float64, 2)
  wvs = zero(Float64)

  # set variables
  ws = StatsBase.rand(settings.distws)

  # allocate the neighbours of i and j
  for u in 1:S
    support.smMP[1, u + 1] = support.smMP[1, u]

    wv1[1] = loglik(likelihood, likelihood.x[S[h]], μ1, τ1)
    wv1[2] = loglik(likelihood, likelihood.x[S[h]], μ2, τ2)

    # wv2[1] = (ws * l(μ1, τ1)) / (ws * l(μ1, τ1) + (1 - ws) * l(μ2, τ2))
    wv2[1] = 1 / (1 + exp(log(1 - ws) + wv1[2] - log(ws) - wv1[1]))
    wv2[2] = 1.0 - wv2[1]

    z = WeightVec(wv2, 1.0)

    pS[h] = sample(1:2, z)

    o = m[pS[h]]
    n[pS[h]] += 1
    m[pS[h]] += (likelihood.x[S[h]] - o) / n[pS[h]]
    d[pS[h]] += (likelihood.x[S[h]] - o) * (likelihood.x[S[h]] - m[pS[h]])

    ll[pS[h]] += wv1[pS[h]]
    ll[3] += loglik(likelihood, likelihood.x[S[h]], μ, τ)

    lq += log(values(z)[pS[h]])
  end

  dist_wm = Beta(splitmerge.para_wm + n[1], splitmerge.para_wm + n[2])

  ratio = exp(split_lr_prior_partition(prior_part, n) +
              split_lr_prior_parameter(prior_para, μ, τ, μ1, τ1, μ2, τ2) +
              split_lr_likelihood(ll) +
              split_lr_operator(splitmerge, dist_wm, ws, u1, u2, lpij, lq) +
              split_ljacobian(τ, τ1, τ2, μ1, μ2, u1, u2))

  if ratio >= 1.0 || rand(Bernoulli(ratio)) == 1
    # move units to their new cluster
    k = findfirst(mcmc_run.empty_clusters)

    mcmc_run.p[ij[2]] = k
    mcmc_run.p[S[pS .== 2]] = k

    # update parameters
    mcmc_run.μ[g] = μ1
    mcmc_run.τ[g] = τ1
    mcmc_run.μ[k] = μ2
    mcmc_run.τ[k] = τ2

    # update statistics
    mcmc_run.s[1, g] = n[1]
    mcmc_run.s[2, g] = m[1]
    mcmc_run.s[3, g] = d[1]
    mcmc_run.s[1, k] = n[2]
    mcmc_run.s[2, k] = m[2]
    mcmc_run.s[3, k] = d[2]

    # set empty cluster
    mcmc_run.empty_clusters[k] = false

    mcmc_run.counter[1, 1] += 1
  end

  mcmc_run.counter[1, 2] += 1

  nothing
end
