# This file is part of Kpax3. License is MIT.

function computeclusteriseqprobs!(x::UInt8,
                                  b::Int,
                                  priorC::AminoAcidPriorCol,
                                  support::MCMCSupport)
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
                                  support::MCMCSupport)
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
                                support::MCMCSupport)
  support.ni[b] += float(x)

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

function updateclusterjweights!(x::UInt8,
                                b::Int,
                                support::MCMCSupport)
  support.nj[b] += float(x)

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

function updateclusteri!(u::Int,
                         data::Matrix{UInt8},
                         support::MCMCSupport)
  support.vi += 1
  support.ui[support.vi] = u

  for b in 1:support.m
    updateclusteriweights!(data[b, u], b, support)
  end

  nothing
end

function updateclusterj!(u::Int,
                         data::Matrix{UInt8},
                         support::MCMCSupport)
  support.vj += 1
  support.uj[support.vj] = u

  for b in 1:support.m
    updateclusterjweights!(data[b, u], b, support)
  end

  nothing
end

function updatelogmargliki!(priorC::AminoAcidPriorCol,
                            support::MCMCSupport)
  for b in 1:support.m
    support.lpi[1, b] = logmarglik(support.ni[b], support.vi, priorC.A[1, b],
                                   priorC.B[1, b])
    support.lpi[2, b] = logmarglik(support.ni[b], support.vi, priorC.A[2, b],
                                   priorC.B[2, b])
    support.lpi[3, b] = logmarglik(support.ni[b], support.vi, priorC.A[3, b],
                                   priorC.B[3, b])
    support.lpi[4, b] = logmarglik(support.ni[b], support.vi, priorC.A[4, b],
                                   priorC.B[4, b])
  end

  nothing
end

function updatelogmargliki!(ni::Vector{Float64},
                            vi::Int,
                            priorC::AminoAcidPriorCol,
                            support::MCMCSupport)
  for b in 1:support.m
    support.lpi[1, b] = logmarglik(ni[b], vi, priorC.A[1, b], priorC.B[1, b])
    support.lpi[2, b] = logmarglik(ni[b], vi, priorC.A[2, b], priorC.B[2, b])
    support.lpi[3, b] = logmarglik(ni[b], vi, priorC.A[3, b], priorC.B[3, b])
    support.lpi[4, b] = logmarglik(ni[b], vi, priorC.A[4, b], priorC.B[4, b])
  end

  nothing
end

function updatelogmarglikj!(priorC::AminoAcidPriorCol,
                            support::MCMCSupport)
  for b in 1:support.m
    support.lpj[1, b] = logmarglik(support.nj[b], support.vj, priorC.A[1, b],
                                   priorC.B[1, b])
    support.lpj[2, b] = logmarglik(support.nj[b], support.vj, priorC.A[2, b],
                                   priorC.B[2, b])
    support.lpj[3, b] = logmarglik(support.nj[b], support.vj, priorC.A[3, b],
                                   priorC.B[3, b])
    support.lpj[4, b] = logmarglik(support.nj[b], support.vj, priorC.A[4, b],
                                   priorC.B[4, b])
  end

  nothing
end

function updatelogmarglikj!(nj::Vector{Float64},
                            vj::Int,
                            priorC::AminoAcidPriorCol,
                            support::MCMCSupport)
  for b in 1:support.m
    support.lpj[1, b] = logmarglik(nj[b], vj, priorC.A[1, b], priorC.B[1, b])
    support.lpj[2, b] = logmarglik(nj[b], vj, priorC.A[2, b], priorC.B[2, b])
    support.lpj[3, b] = logmarglik(nj[b], vj, priorC.A[3, b], priorC.B[3, b])
    support.lpj[4, b] = logmarglik(nj[b], vj, priorC.A[4, b], priorC.B[4, b])
  end

  nothing
end

function initclusteriweights!(x::UInt8,
                              b::Int,
                              k::Int,
                              priorC::AminoAcidPriorCol,
                              support::MCMCSupport)
  support.ni[b] = float(x)

  support.wi.w[1, b] = priorC.logγ[1] +
                       logmarglik(x, 1, priorC.A[1, b], priorC.B[1, b])
  support.wi.w[2, b] = priorC.logγ[2] +
                       logmarglik(x, 1, priorC.A[2, b], priorC.B[2, b])
  support.wi.w[3, b] = priorC.logγ[3] + priorC.logω[k][1] +
                       logmarglik(x, 1, priorC.A[3, b], priorC.B[3, b])
  support.wi.w[4, b] = priorC.logγ[3] + priorC.logω[k][2] +
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
                              k::Int,
                              priorC::AminoAcidPriorCol,
                              support::MCMCSupport)
  support.nj[b] = float(x)

  support.wj.w[1, b] = priorC.logγ[1] +
                       logmarglik(x, 1, priorC.A[1, b], priorC.B[1, b])
  support.wj.w[2, b] = priorC.logγ[2] +
                       logmarglik(x, 1, priorC.A[2, b], priorC.B[2, b])
  support.wj.w[3, b] = priorC.logγ[3] + priorC.logω[k][1] +
                       logmarglik(x, 1, priorC.A[3, b], priorC.B[3, b])
  support.wj.w[4, b] = priorC.logγ[3] + priorC.logω[k][2] +
                       logmarglik(x, 1, priorC.A[4, b], priorC.B[4, b])

  M = max(support.wj.w[1, b], support.wj.w[2, b],
          support.wj.w[3, b], support.wj.w[4, b])

  support.wj.c[b] = M + log(exp(support.wj.w[1, b] - M) +
                            exp(support.wj.w[2, b] - M) +
                            exp(support.wj.w[3, b] - M) +
                            exp(support.wj.w[4, b] - M))

  nothing
end

function initsupportsplit!(ij::Vector{Int},
                           k::Int,
                           data::Matrix{UInt8},
                           priorC::AminoAcidPriorCol,
                           settings::KSettings,
                           support::MCMCSupport)
  len = min(support.n, k + settings.maxclust - 1)

  resizelogω!(priorC, len)
  resizesupport!(support, len)

  support.vi = 1
  support.vj = 1

  support.ui[1] = ij[1]
  support.uj[1] = ij[2]

  for b in 1:support.m
    initclusteriweights!(data[b, ij[1]], b, k, priorC, support)
    initclusterjweights!(data[b, ij[2]], b, k, priorC, support)
  end

  nothing
end

function initsupportmerge!(ij::Vector{Int},
                           k::Int,
                           data::Matrix{UInt8},
                           priorC::AminoAcidPriorCol,
                           support::MCMCSupport)
  support.vi = 1
  support.vj = 1

  support.ui[1] = ij[1]
  support.uj[1] = ij[2]

  for b in 1:support.m
    initclusteriweights!(data[b, ij[1]], b, k, priorC, support)
    initclusterjweights!(data[b, ij[2]], b, k, priorC, support)
  end

  nothing
end

function initsupportbrwsplit!(k::Int,
                              i::Int,
                              hi::Int,
                              data::Matrix{UInt8},
                              priorC::AminoAcidPriorCol,
                              settings::KSettings,
                              support::MCMCSupport,
                              state::AminoAcidState)
  len = min(support.n, k + settings.maxclust - 1)

  resizelogω!(priorC, len)
  resizesupport!(support, len)

  support.vi = state.v[hi] - 1
  support.vj = 1

  for b in 1:support.m
    support.nj[b] = float(data[b, i])
    support.ni[b] = state.n1s[hi, b] - support.nj[b]
  end

  idx = 0
  for a in 1:state.v[hi]
    if state.unit[hi][a] != i
      idx += 1
      support.ui[idx] = state.unit[hi][a]
    end
  end

  nothing
end

function initsupportbrwmerge!(i::Int,
                              hj::Int,
                              data::Matrix{UInt8},
                              support::MCMCSupport,
                              state::AminoAcidState)
  support.vj = state.v[hj] + 1

  for b in 1:support.m
    support.nj[b] = state.n1s[hj, b] + float(data[b, i])
  end

  nothing
end

function initsupportbrwmove!(i::Int,
                             hi::Int,
                             hj::Int,
                             data::Matrix{UInt8},
                             support::MCMCSupport,
                             state::AminoAcidState)
  support.vi = state.v[hi] - 1
  support.vj = state.v[hj] + 1

  for b in 1:support.m
    support.ni[b] = state.n1s[hi, b] - float(data[b, i])
    support.nj[b] = state.n1s[hj, b] + float(data[b, i])
  end

  idx = 0
  for a in 1:state.v[hi]
    if state.unit[hi][a] != i
      idx += 1
      support.ui[idx] = state.unit[hi][a]
    end
  end

  nothing
end
