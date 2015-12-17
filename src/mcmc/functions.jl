# This file is part of Kpax3. License is MIT.

function computeclusteriseqprobs!(x::UInt8,
                                  b::Int,
                                  priorC::AminoAcidPriorCol,
                                  support::KSupport)
  support.wi.z[1, b] = support.wi.w[1, b] +
                       logcondmarglik(x, support.ni[b], support.vi,
                                      priorC.A[1, b], priorC.B[1, b])
  support.wi.z[2, b] = support.wi.w[2, b] +
                       logcondmarglik(x, support.ni[b], support.vi,
                                      priorC.A[2, b], priorC.B[2, b])
  support.wi.z[3, b] = support.wi.w[3, b] +
                       logcondmarglik(x, support.ni[b], support.vi,
                                      priorC.A[3, b], priorC.B[3, b])
  support.wi.z[4, b] = support.wi.w[4, b] +
                       logcondmarglik(x, support.ni[b], support.vi,
                                      priorC.A[4, b], priorC.B[4, b])

  support.tmp[1] = support.wi.z[1, b] - support.wi.c[b]
  support.tmp[2] = support.wi.z[2, b] - support.wi.c[b]
  support.tmp[3] = support.wi.z[3, b] - support.wi.c[b]
  support.tmp[4] = support.wi.z[4, b] - support.wi.c[b]

  M = max(support.tmp[1], support.tmp[2], support.tmp[3], support.tmp[4])

  M + log(exp(support.tmp[1] - M) + exp(support.tmp[2] - M) +
          exp(support.tmp[3] - M) + exp(support.tmp[4] - M))
end

function computeclusterjseqprobs!(x::UInt8,
                                  b::Int,
                                  priorC::AminoAcidPriorCol,
                                  support::KSupport)
  support.wj.z[1, b] = support.wj.w[1, b] +
                       logcondmarglik(x, support.nj[b], support.vj,
                                      priorC.A[1, b], priorC.B[1, b])
  support.wj.z[2, b] = support.wj.w[2, b] +
                       logcondmarglik(x, support.nj[b], support.vj,
                                      priorC.A[2, b], priorC.B[2, b])
  support.wj.z[3, b] = support.wj.w[3, b] +
                       logcondmarglik(x, support.nj[b], support.vj,
                                      priorC.A[3, b], priorC.B[3, b])
  support.wj.z[4, b] = support.wj.w[4, b] +
                       logcondmarglik(x, support.nj[b], support.vj,
                                      priorC.A[4, b], priorC.B[4, b])

  support.tmp[1] = support.wj.z[1, b] - support.wj.c[b]
  support.tmp[2] = support.wj.z[2, b] - support.wj.c[b]
  support.tmp[3] = support.wj.z[3, b] - support.wj.c[b]
  support.tmp[4] = support.wj.z[4, b] - support.wj.c[b]

  M = max(support.tmp[1], support.tmp[2], support.tmp[3], support.tmp[4])

  M + log(exp(support.tmp[1] - M) + exp(support.tmp[2] - M) +
          exp(support.tmp[3] - M) + exp(support.tmp[4] - M))
end

function updateclusteriweights!(x::UInt8,
                                b::Int,
                                support::KSupport)
  support.ni[b] += Float64(x)

  support.wi.w[1, b] = support.wi.z[1, b]
  support.wi.w[2, b] = support.wi.z[2, b]
  support.wi.w[3, b] = support.wi.z[3, b]
  support.wi.w[4, b] = support.wi.z[4, b]

  M = max(support.wi.w[1, b], support.wi.w[2, b], support.wi.w[3, b],
          support.wi.w[4, b])

  support.wi.c[b] = M + log(exp(support.wi.w[1, b] - M) +
                            exp(support.wi.w[2, b] - M) +
                            exp(support.wi.w[3, b] - M) +
                            exp(support.wi.w[4, b] - M))

  nothing
end

function updateclusteri!(u::Int,
                         data::Matrix{UInt8},
                         support::KSupport)
  support.vi += 1
  support.ui[support.vi] = u

  for b in 1:size(data, 1)
    updateclusteriweights!(data[b, u], b, support)
  end

  nothing
end

function updateclusterjweights!(x::UInt8,
                                b::Int,
                                support::KSupport)
  support.nj[b] += Float64(x)

  support.wj.w[1, b] = support.wj.z[1, b]
  support.wj.w[2, b] = support.wj.z[2, b]
  support.wj.w[3, b] = support.wj.z[3, b]
  support.wj.w[4, b] = support.wj.z[4, b]

  M = max(support.wj.w[1, b], support.wj.w[2, b], support.wj.w[3, b],
          support.wj.w[4, b])

  support.wj.c[b] = M + log(exp(support.wj.w[1, b] - M) +
                            exp(support.wj.w[2, b] - M) +
                            exp(support.wj.w[3, b] - M) +
                            exp(support.wj.w[4, b] - M))

  nothing
end

function updateclusterj!(u::Int,
                         data::Matrix{UInt8},
                         support::KSupport)
  support.vj += 1
  support.uj[support.vj] = u

  for b in 1:size(data, 1)
    updateclusterjweights!(data[b, u], b, support)
  end

  nothing
end

function initclusteriweights!(x::UInt8,
                              b::Int,
                              logω::Vector{Float64},
                              priorC::AminoAcidPriorCol,
                              support::KSupport)
  support.ni[b] = Float64(x)

  support.wi.w[1, b] = priorC.logγ[1] + logω[1] +
                       logmarglik(x, 1, priorC.A[1, b], priorC.B[1, b])
  support.wi.w[2, b] = priorC.logγ[2] + logω[2] +
                       logmarglik(x, 1, priorC.A[2, b], priorC.B[2, b])
  support.wi.w[3, b] = priorC.logγ[3] + logω[3] +
                       logmarglik(x, 1, priorC.A[3, b], priorC.B[3, b])
  support.wi.w[4, b] = priorC.logγ[4] + logω[4] +
                       logmarglik(x, 1, priorC.A[4, b], priorC.B[4, b])

  M = max(support.wi.w[1, b], support.wi.w[2, b],
          support.wi.w[3, b], support.wi.w[4, b])

  support.wi.c[b] = M + log(exp(support.wi.w[1, b] - M) +
                            exp(support.wi.w[2, b] - M) +
                            exp(support.wi.w[3, b] - M) +
                            exp(support.wi.w[4, b] - M))

  nothing
end

function initclusterjweights!(x::UInt8,
                              b::Int,
                              logω::Vector{Float64},
                              priorC::AminoAcidPriorCol,
                              support::KSupport)
  support.nj[b] = Float64(x)

  support.wj.w[1, b] = priorC.logγ[1] + logω[1] +
                       logmarglik(x, 1, priorC.A[1, b], priorC.B[1, b])
  support.wj.w[2, b] = priorC.logγ[2] + logω[2] +
                       logmarglik(x, 1, priorC.A[2, b], priorC.B[2, b])
  support.wj.w[3, b] = priorC.logγ[3] + logω[3] +
                       logmarglik(x, 1, priorC.A[3, b], priorC.B[3, b])
  support.wj.w[4, b] = priorC.logγ[4] + logω[4] +
                       logmarglik(x, 1, priorC.A[4, b], priorC.B[4, b])

  M = max(support.wj.w[1, b], support.wj.w[2, b],
          support.wj.w[3, b], support.wj.w[4, b])

  support.wj.c[b] = M + log(exp(support.wj.w[1, b] - M) +
                            exp(support.wj.w[2, b] - M) +
                            exp(support.wj.w[3, b] - M) +
                            exp(support.wj.w[4, b] - M))

  nothing
end

function initsupport!(ij::Vector{Int},
                      S::Int,
                      k::Int,
                      data::Matrix{UInt8},
                      logω::Vector{Float64},
                      priorC::AminoAcidPriorCol,
                      support::KSupport)
  support.vi = 1
  support.vj = 1

  if length(support.ui) <= S
    support.ui = zeros(Int, S + 1)
    support.uj = zeros(Int, S + 1)
  end

  support.ui[1] = ij[1]
  support.uj[1] = ij[2]

  if size(support.C, 1) < k
    support.C = zeros(UInt8, k, size(data, 1))
  end

  for b in 1:size(data, 1)
    initclusteriweights!(data[b, ij[1]], b, logω, priorC, support)
    initclusterjweights!(data[b, ij[2]], b, logω, priorC, support)
  end

  nothing
end
