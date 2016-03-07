# This file is part of Kpax3. License is MIT.

function samplenoise!(C::Matrix{UInt8},
                      logp::Vector{Float64},
                      cl,
                      b::Int,
                      x::Float64,
                      logγ::Vector{Float64})
  logp[1] += logγ[1]
  logp[2] += x

  for g in cl
    C[g, b] = 0x01
  end

  nothing
end

function sampleweaksignal!(C::Matrix{UInt8},
                           logp::Vector{Float64},
                           cl,
                           b::Int,
                           x::Float64,
                           logγ::Vector{Float64})
  logp[1] += logγ[2]
  logp[2] += x

  for g in cl
    C[g, b] = 0x02
  end

  nothing
end

function samplestrongsignal!(C::Matrix{UInt8},
                             logp::Vector{Float64},
                             cl,
                             b::Int,
                             x::Float64,
                             logγ::Vector{Float64},
                             logω::Vector{Float64},
                             lω::Matrix{Float64})
  logp[1] += logγ[3]
  logp[2] += x

  for l in 1:length(cl)
    if rand() <= exp(lω[1, l])
      C[cl[l], b] = 0x03

      logp[1] += logω[3]
      logp[2] += lω[1, l]
    else
      C[cl[l], b] = 0x04

      logp[1] += logω[4]
      logp[2] += lω[2, l]
    end
  end

  nothing
end

function computeclusterlogprobs!(logq::Matrix{Float64},
                                 lγ::Vector{Float64},
                                 lω::Matrix{Float64},
                                 b::Int,
                                 l::Int,
                                 y::Real,
                                 n::Real,
                                 logω::Vector{Float64},
                                 priorC::AminoAcidPriorCol)
  logq[1, l] = logmarglik(y, n, priorC.A[1, b], priorC.B[1, b])
  logq[2, l] = logmarglik(y, n, priorC.A[2, b], priorC.B[2, b])
  logq[3, l] = logω[3] + logmarglik(y, n, priorC.A[3, b], priorC.B[3, b])
  logq[4, l] = logω[4] + logmarglik(y, n, priorC.A[4, b], priorC.B[4, b])

  lγ[1] += logq[1, l]
  lγ[2] += logq[2, l]

  if logq[3, l] > logq[4, l]
    tmp = log1p(exp(logq[4, l] - logq[3, l]))

    lγ[3] += logq[3, l] + tmp

    lω[1, l] = -tmp
    lω[2, l] = logq[4, l] - logq[3, l] - tmp
  else
    tmp = log1p(exp(logq[3, l] - logq[4, l]))

    lγ[3] += logq[4, l] + tmp

    lω[1, l] = logq[3, l] - logq[4, l] - tmp
    lω[2, l] = -tmp
  end

  nothing
end

function rpostpartitioncols!(C::Matrix{UInt8},
                             cl::Vector{Int},
                             v::Vector{Int},
                             n1s::Matrix{Float64},
                             priorC::AminoAcidPriorCol)
  logp = zeros(Float64, 2)

  logq = zeros(Float64, 4, length(cl))

  lγ = zeros(Float64, 3)
  lω = zeros(Float64, 2, length(cl))

  g = 0
  p = 0.0
  tmp = 0.0

  for b in 1:size(C, 2)
    lγ[1] = priorC.logγ[1]
    lγ[2] = priorC.logγ[2]
    lγ[3] = priorC.logγ[3]

    for l in 1:length(cl)
      g = cl[l]
      computeclusterlogprobs!(logq, lγ, lω, b, l, n1s[g, b], v[g], priorC.logω,
                              priorC)
    end

    if (lγ[1] >= lγ[2]) && (lγ[1] >= lγ[3])
      tmp = lγ[1] + log1p(exp(lγ[2] - lγ[1]) + exp(lγ[3] - lγ[1]))
    elseif (lγ[2] >= lγ[1]) && (lγ[2] >= lγ[3])
      tmp = lγ[2] + log1p(exp(lγ[1] - lγ[2]) + exp(lγ[3] - lγ[2]))
    else
      tmp = lγ[3] + log1p(exp(lγ[1] - lγ[3]) + exp(lγ[2] - lγ[3]))
    end

    lγ[1] -= tmp
    lγ[2] -= tmp
    lγ[3] -= tmp

    p = rand()

    if p <= exp(lγ[1])
      samplenoise!(C, logp, cl, b, lγ[1], priorC.logγ)
    elseif p <= exp(lγ[1]) + exp(lγ[2])
      sampleweaksignal!(C, logp, cl, b, lγ[2], priorC.logγ)
    else
      samplestrongsignal!(C, logp, cl, b, lγ[3], priorC.logγ, priorC.logω, lω)
    end
  end

  logp
end

function simcbrw!(k::Int,
                  hi::Int,
                  hj::Int,
                  priorC::PriorColPartition,
                  support::KSupport,
                  mcmcobj::AminoAcidMCMC)
  fill!(support.logpC, 0.0)

  logq = zeros(Float64, 4, k)

  lγ = zeros(Float64, 3)
  lω = zeros(Float64, 2, k)

  p = 0.0
  tmp = 0.0

  for b in 1:size(support.C, 2)
    lγ[1] = priorC.logγ[1]
    lγ[2] = priorC.logγ[2]
    lγ[3] = priorC.logγ[3]

    l = 1
    for g in mcmcobj.cl
      if g != hi
        if g != hj
          # this cluster is not being affected by the move
          computeclusterlogprobs!(logq, lγ, lω, b, l, mcmcobj.n1s[g, b],
                                  mcmcobj.v[g], support.logω, priorC)
        else
          # we are moving i into an existing cluster
          computeclusterlogprobs!(logq, lγ, lω, b, l,
                                  mcmcobj.n1s[g, b] + support.ni[b],
                                  mcmcobj.v[g] + 1, support.logω, priorC)
        end

        l += 1
      else
        if g == hj
          # move singleton i into the same cluster
          computeclusterlogprobs!(logq, lγ, lω, b, l, support.ni[b], 1,
                                  support.logω, priorC)
          l += 1
        elseif mcmcobj.v[g] > 1
          # remove i from this cluster
          computeclusterlogprobs!(logq, lγ, lω, b, l,
                                  mcmcobj.n1s[g, b] - support.ni[b],
                                  mcmcobj.v[g] - 1, support.logω, priorC)
          l += 1
        end
      end
    end

    if length(mcmcobj.cl) < k
      # move i into a new cluster of its own
      computeclusterlogprobs!(logq, lγ, lω, b, l, support.ni[b], 1,
                              support.logω, priorC)
    end

    if (lγ[1] >= lγ[2]) && (lγ[1] >= lγ[3])
      tmp = lγ[1] + log1p(exp(lγ[2] - lγ[1]) + exp(lγ[3] - lγ[1]))
    elseif (lγ[2] >= lγ[1]) && (lγ[2] >= lγ[3])
      tmp = lγ[2] + log1p(exp(lγ[1] - lγ[2]) + exp(lγ[3] - lγ[2]))
    else
      tmp = lγ[3] + log1p(exp(lγ[1] - lγ[3]) + exp(lγ[2] - lγ[3]))
    end

    lγ[1] -= tmp
    lγ[2] -= tmp
    lγ[3] -= tmp

    p = rand()

    if p <= exp(lγ[1])
      samplenoise!(support.C, support.logpC, 1:k, b, lγ[1], priorC.logγ)
    elseif p <= exp(lγ[1]) + exp(lγ[2])
      sampleweaksignal!(support.C, support.logpC, 1:k, b, lγ[2], priorC.logγ)
    else
      samplestrongsignal!(support.C, support.logpC, 1:k, b, lγ[3], priorC.logγ,
                          support.logω, lω)
    end
  end

  nothing
end

function simcmerge!(k::Int,
                    hi::Int,
                    hj::Int,
                    vi::Int,
                    ni::Vector{Float64},
                    priorC::PriorColPartition,
                    support::KSupport,
                    mcmcobj::AminoAcidMCMC)
  fill!(support.logpC, 0.0)

  logq = zeros(Float64, 4, k)

  lγ = zeros(Float64, 3)
  lω = zeros(Float64, 2, k)

  p = 0.0
  tmp = 0.0

  for b in 1:size(support.C, 2)
    lγ[1] = priorC.logγ[1]
    lγ[2] = priorC.logγ[2]
    lγ[3] = priorC.logγ[3]

    l = 1
    for g in mcmcobj.cl
      if g != hj
        if g != hi
          computeclusterlogprobs!(logq, lγ, lω, b, l, mcmcobj.n1s[g, b],
                                  mcmcobj.v[g], support.logω, priorC)
        else
          computeclusterlogprobs!(logq, lγ, lω, b, l, ni[b], vi, support.logω,
                                  priorC)
        end

        l += 1
      end
    end

    if (lγ[1] >= lγ[2]) && (lγ[1] >= lγ[3])
      tmp = lγ[1] + log1p(exp(lγ[2] - lγ[1]) + exp(lγ[3] - lγ[1]))
    elseif (lγ[2] >= lγ[1]) && (lγ[2] >= lγ[3])
      tmp = lγ[2] + log1p(exp(lγ[1] - lγ[2]) + exp(lγ[3] - lγ[2]))
    else
      tmp = lγ[3] + log1p(exp(lγ[1] - lγ[3]) + exp(lγ[2] - lγ[3]))
    end

    lγ[1] -= tmp
    lγ[2] -= tmp
    lγ[3] -= tmp

    p = rand()

    if p <= exp(lγ[1])
      samplenoise!(support.C, support.logpC, 1:k, b, lγ[1], priorC.logγ)
    elseif p <= exp(lγ[1]) + exp(lγ[2])
      sampleweaksignal!(support.C, support.logpC, 1:k, b, lγ[2], priorC.logγ)
    else
      samplestrongsignal!(support.C, support.logpC, 1:k, b, lγ[3], priorC.logγ,
                          support.logω, lω)
    end
  end

  nothing
end

function simcsplit!(k::Int,
                    hi::Int,
                    priorC::PriorColPartition,
                    support::KSupport,
                    mcmcobj::AminoAcidMCMC)
  fill!(support.logpC, 0.0)

  logq = zeros(Float64, 4, k)

  lγ = zeros(Float64, 3)
  lω = zeros(Float64, 2, k)

  p = 0.0
  tmp = 0.0

  for b in 1:size(support.C, 2)
    lγ[1] = priorC.logγ[1]
    lγ[2] = priorC.logγ[2]
    lγ[3] = priorC.logγ[3]

    l = 1
    for g in mcmcobj.cl
      if g != hi
        computeclusterlogprobs!(logq, lγ, lω, b, l, mcmcobj.n1s[g, b],
                                mcmcobj.v[g], support.logω, priorC)
      else
        computeclusterlogprobs!(logq, lγ, lω, b, l, support.ni[b], support.vi,
                                support.logω, priorC)
      end

      l += 1
    end

    computeclusterlogprobs!(logq, lγ, lω, b, l, support.nj[b], support.vj,
                            support.logω, priorC)

    if (lγ[1] >= lγ[2]) && (lγ[1] >= lγ[3])
      tmp = lγ[1] + log1p(exp(lγ[2] - lγ[1]) + exp(lγ[3] - lγ[1]))
    elseif (lγ[2] >= lγ[1]) && (lγ[2] >= lγ[3])
      tmp = lγ[2] + log1p(exp(lγ[1] - lγ[2]) + exp(lγ[3] - lγ[2]))
    else
      tmp = lγ[3] + log1p(exp(lγ[1] - lγ[3]) + exp(lγ[2] - lγ[3]))
    end

    lγ[1] -= tmp
    lγ[2] -= tmp
    lγ[3] -= tmp

    p = rand()

    if p <= exp(lγ[1])
      samplenoise!(support.C, support.logpC, 1:k, b, lγ[1], priorC.logγ)
    elseif p <= exp(lγ[1]) + exp(lγ[2])
      sampleweaksignal!(support.C, support.logpC, 1:k, b, lγ[2], priorC.logγ)
    else
      samplestrongsignal!(support.C, support.logpC, 1:k, b, lγ[3], priorC.logγ,
                          support.logω, lω)
    end
  end

  nothing
end
