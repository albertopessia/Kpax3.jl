# This file is part of Kpax3. License is MIT.

function initializeSupport!(support::KSupport,
                            i::Int,
                            xi::Array{UInt8, 1},
                            j::Int,
                            xj::Array{UInt8, 1},
                            S::Int,
                            logγ::Array{Float64, 1},
                            logω::Array{Float64, 1},
                            A::Array{Float64, 2},
                            B::Array{Float64, 2})
  M = zero(Float64)

  support.gi.v = one(Int)
  if length(support.gi.unit) < S + 1
    support.gi.unit = zeros(Int, S + 1)
  end
  support.gi.unit[1] = i

  support.gj.v = one(Int)
  if length(support.gj.unit) < S + 1
    support.gj.unit = zeros(Int, S + 1)
  end
  support.gj.unit[1] = j

  for b in 1:length(xi)
    support.gi.n1s[b] = Float64(xi[b])

    support.wi.w[1, b] = logγ[1] + logω[1] +
                         logmarglik(xi[b], 1, A[1, b], B[1, b])
    support.wi.w[2, b] = logγ[2] + logω[2] +
                         logmarglik(xi[b], 1, A[2, b], B[2, b])
    support.wi.w[3, b] = logγ[3] + logω[3] +
                         logmarglik(xi[b], 1, A[3, b], B[3, b])
    support.wi.w[4, b] = logγ[4] + logω[4] +
                         logmarglik(xi[b], 1, A[4, b], B[4, b])

    M = max(support.wi.w[1, b], support.wi.w[2, b],
            support.wi.w[3, b], support.wi.w[4, b])

    support.wi.c[b] = M + log(exp(support.wi.w[1, b] - M) +
                              exp(support.wi.w[2, b] - M) +
                              exp(support.wi.w[3, b] - M) +
                              exp(support.wi.w[4, b] - M))

    support.gj.n1s[b] = Float64(xj[b])

    support.wj.w[1, b] = logγ[1] + logω[1] +
                         logmarglik(xj[b], 1, A[1, b], B[1, b])
    support.wj.w[2, b] = logγ[2] + logω[2] +
                         logmarglik(xj[b], 1, A[2, b], B[2, b])
    support.wj.w[3, b] = logγ[3] + logω[3] +
                         logmarglik(xj[b], 1, A[3, b], B[3, b])
    support.wj.w[4, b] = logγ[4] + logω[4] +
                         logmarglik(xj[b], 1, A[4, b], B[4, b])

    M = max(support.wj.w[1, b], support.wj.w[2, b],
            support.wj.w[3, b], support.wj.w[4, b])

    support.wj.c[b] = M + log(exp(support.wj.w[1, b] - M) +
                              exp(support.wj.w[2, b] - M) +
                              exp(support.wj.w[3, b] - M) +
                              exp(support.wj.w[4, b] - M))
  end

  nothing
end

function updateSupport!(h::KCluster,
                        w::KWeight,
                        u::Int,
                        x::Array{UInt8, 1})
  M = zero(Float64)

  h.v += 1
  h.unit[h.v] = u

  for b in 1:length(x)
    h.n1s[b] += Float64(x[b])

    w.w[1, b] = w.z[1, b]
    w.w[2, b] = w.z[2, b]
    w.w[3, b] = w.z[3, b]
    w.w[4, b] = w.z[4, b]

    M = max(w.w[1, b], w.w[2, b], w.w[3, b], w.w[4, b])

    w.c[b] = M + log(exp(w.w[1, b] - M) + exp(w.w[2, b] - M) +
                     exp(w.w[3, b] - M) + exp(w.w[4, b] - M))
  end

  nothing
end
