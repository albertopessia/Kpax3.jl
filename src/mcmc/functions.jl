# This file is part of Kpax3. License is MIT.

function initializeSupport!(support::KSupport,
                            i::Int,
                            j::Int,
                            S::Int,
                            k::Int,
                            data::Matrix{UInt8},
                            logγ::Vector{Float64},
                            logω::Vector{Float64},
                            A::Matrix{Float64},
                            B::Matrix{Float64})
  M = 0.0

  support.vi = 1
  if length(support.ui) <= S
    support.ui = zeros(Int, S + 1)
  end
  support.ui[1] = i

  support.vj = 1
  if length(support.uj) <= S
    support.uj = zeros(Int, S + 1)
  end
  support.uj[1] = j

  if length(support.filledcluster) >= k
    fill!(support.C, 0x00)
    fill!(support.filledcluster, false)
    fill!(support.v, 0)
    fill!(support.n1s, 0.0)
  else
    support.C = zeros(UInt8, k, size(data, 1))
    support.filledcluster = falses(k)
    support.v = zeros(Int, k)
    support.n1s = zeros(Float64, k, size(data, 1))
  end

  for b in 1:size(data, 1)
    support.ni[b] = Float64(data[b, i])

    support.wi.w[1, b] = logγ[1] + logω[1] +
                         logmarglik(data[b, i], 1, A[1, b], B[1, b])
    support.wi.w[2, b] = logγ[2] + logω[2] +
                         logmarglik(data[b, i], 1, A[2, b], B[2, b])
    support.wi.w[3, b] = logγ[3] + logω[3] +
                         logmarglik(data[b, i], 1, A[3, b], B[3, b])
    support.wi.w[4, b] = logγ[4] + logω[4] +
                         logmarglik(data[b, i], 1, A[4, b], B[4, b])

    M = max(support.wi.w[1, b], support.wi.w[2, b],
            support.wi.w[3, b], support.wi.w[4, b])

    support.wi.c[b] = M + log(exp(support.wi.w[1, b] - M) +
                              exp(support.wi.w[2, b] - M) +
                              exp(support.wi.w[3, b] - M) +
                              exp(support.wi.w[4, b] - M))

    support.nj[b] = Float64(data[b, j])

    support.wj.w[1, b] = logγ[1] + logω[1] +
                         logmarglik(data[b, j], 1, A[1, b], B[1, b])
    support.wj.w[2, b] = logγ[2] + logω[2] +
                         logmarglik(data[b, j], 1, A[2, b], B[2, b])
    support.wj.w[3, b] = logγ[3] + logω[3] +
                         logmarglik(data[b, j], 1, A[3, b], B[3, b])
    support.wj.w[4, b] = logγ[4] + logω[4] +
                         logmarglik(data[b, j], 1, A[4, b], B[4, b])

    M = max(support.wj.w[1, b], support.wj.w[2, b],
            support.wj.w[3, b], support.wj.w[4, b])

    support.wj.c[b] = M + log(exp(support.wj.w[1, b] - M) +
                              exp(support.wj.w[2, b] - M) +
                              exp(support.wj.w[3, b] - M) +
                              exp(support.wj.w[4, b] - M))
  end

  nothing
end

function updateClusteri!(support::KSupport,
                         u::Int,
                         data::Matrix{UInt8})
  M = 0.0

  support.vi += 1
  support.ui[support.vi] = u

  for b in 1:size(data, 1)
    support.ni[b] += Float64(data[b, u])

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
  end

  nothing
end

function updateClusterj!(support::KSupport,
                         u::Int,
                         data::Matrix{UInt8})
  M = 0.0

  support.vj += 1
  support.uj[support.vj] = u

  for b in 1:size(data, 1)
    support.nj[b] += Float64(data[b, u])

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
  end

  nothing
end
