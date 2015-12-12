# This file is part of Kpax3. License is MIT.

function rpostpartitioncols!(C::Matrix{UInt8},
                             cl::Vector{Int},
                             v::Vector{Int},
                             n1s::Matrix{Float64},
                             logγ::Vector{Float64},
                             logω::Vector{Float64},
                             A::Matrix{Float64},
                             B::Matrix{Float64})
  logpr = 0.0
  logpo = 0.0

  logq = zeros(Float64, 4, length(cl))

  lγ = zeros(Float64, 3)
  lω = zeros(Float64, 2, length(cl))

  g = 0
  M = 0.0
  u = 0.0
  tmp = 0.0

  for b in 1:size(C, 2)
    lγ[1] = logγ[1]
    lγ[2] = logγ[2]
    lγ[3] = logγ[3]

    for l in 1:length(cl)
      g = cl[l]

      logq[1, l] = logmarglik(n1s[g, b], v[g], A[1, b], B[1, b])
      logq[2, l] = logmarglik(n1s[g, b], v[g], A[2, b], B[2, b])
      logq[3, l] = logω[3] + logmarglik(n1s[g, b], v[g], A[3, b], B[3, b])
      logq[4, l] = logω[4] + logmarglik(n1s[g, b], v[g], A[4, b], B[4, b])

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
    end

    if (lγ[1] >= lγ[2]) && (lγ[1] >= lγ[3])
      tmp = lγ[1] + log1p(exp(lγ[2] - lγ[1]) + exp(lγ[3] - lγ[1]))
    elseif (lγ[2] >= lγ[1]) && (lγ[2] >= lγ[3])
      tmp = lγ[2] + log1p(exp(lγ[1] - lγ[2]) + exp(lγ[3] - lγ[2]))
    else
      tmp = lγ[3] + log1p(exp(lγ[1] - lγ[3]) + exp(lγ[2] - lγ[3]))
    end

    u = rand()

    if u <= exp(lγ[1] - tmp)
      for g in cl
        C[g, b] = 0x01
      end
      logpr += logγ[1] + logω[1]
      logpo += lγ[1] - tmp
    elseif u <= exp(lγ[1] - tmp) + exp(lγ[2] - tmp)
      for g in cl
        C[g, b] = 0x02
      end
      logpr += logγ[2] + logω[2]
      logpo += lγ[2] - tmp
    else
      for l in 1:length(cl)
        if rand() <= exp(lω[1, l])
          C[cl[l], b] = 0x03
          logpr += logω[3]
          logpo += lω[1, l]
        else
          C[cl[l], b] = 0x04
          logpr += logω[4]
          logpo += lω[2, l]
        end
      end
      logpr += logγ[3] + logω[3]
      logpo += lγ[3] - tmp
    end
  end

  (logpr, logpo)
end

function rpostpartitioncols!(C::Matrix{UInt8},
                             k::Int,
                             hi::Int,
                             vi::Int,
                             ni::Vector{Float64},
                             vj::Int,
                             nj::Vector{Float64},
                             cl::Vector{Int},
                             v::Vector{Int},
                             n1s::Matrix{Float64},
                             logγ::Vector{Float64},
                             logω::Vector{Float64},
                             A::Matrix{Float64},
                             B::Matrix{Float64})
  logpr = 0.0
  logpo = 0.0

  logq = zeros(Float64, 4, k)

  lγ = zeros(Float64, 3)
  lω = zeros(Float64, 2, k)

  g = 0
  M = 0.0
  u = 0.0
  tmp = 0.0

  for b in 1:size(C, 2)
    lγ[1] = logγ[1]
    lγ[2] = logγ[2]
    lγ[3] = logγ[3]

    for l in 1:(k - 1)
      g = cl[l]

      if g != hi
        logq[1, l] = logmarglik(n1s[g, b], v[g], A[1, b], B[1, b])
        logq[2, l] = logmarglik(n1s[g, b], v[g], A[2, b], B[2, b])
        logq[3, l] = logω[3] + logmarglik(n1s[g, b], v[g], A[3, b], B[3, b])
        logq[4, l] = logω[4] + logmarglik(n1s[g, b], v[g], A[4, b], B[4, b])

        lγ[1] += logq[1, l]
        lγ[2] += logq[2, l]

        if logq[3, l] > logq[4, l]
          tmp = log1p(exp(logq[4, l] - logq[3, l]))

          lγ[3] += logq[3, l] + tmp

          lω[1, l] = exp(-tmp)
          lω[2, l] = exp(logq[4, l] - logq[3, l] - tmp)
        else
          tmp = log1p(exp(logq[3, l] - logq[4, l]))

          lγ[3] += logq[4, l] + tmp

          lω[1, l] = exp(logq[3, l] - logq[4, l] - tmp)
          lω[2, l] = exp(-tmp)
        end
      else
        logq[1, l] = logmarglik(ni[b], vi, A[1, b], B[1, b])
        logq[2, l] = logmarglik(ni[b], vi, A[2, b], B[2, b])
        logq[3, l] = logω[3] + logmarglik(ni[b], vi, A[3, b], B[3, b])
        logq[4, l] = logω[4] + logmarglik(ni[b], vi, A[4, b], B[4, b])

        lγ[1] += logq[1, l]
        lγ[2] += logq[2, l]

        if logq[3, l] > logq[4, l]
          tmp = log1p(exp(logq[4, l] - logq[3, l]))

          lγ[3] += logq[3, l] + tmp

          lω[1, l] = exp(-tmp)
          lω[2, l] = exp(logq[4, l] - logq[3, l] - tmp)
        else
          tmp = log1p(exp(logq[3, l] - logq[4, l]))

          lγ[3] += logq[4, l] + tmp

          lω[1, l] = exp(logq[3, l] - logq[4, l] - tmp)
          lω[2, l] = exp(-tmp)
        end
      end
    end

    logq[1, k] = logmarglik(nj[b], vj, A[1, b], B[1, b])
    logq[2, k] = logmarglik(nj[b], vj, A[2, b], B[2, b])
    logq[3, k] = logω[3] + logmarglik(nj[b], vj, A[3, b], B[3, b])
    logq[4, k] = logω[4] + logmarglik(nj[b], vj, A[4, b], B[4, b])

    lγ[1] += logq[1, k]
    lγ[2] += logq[2, k]

    if logq[3, k] > logq[4, k]
      tmp = log1p(exp(logq[4, k] - logq[3, k]))

      lγ[3] += logq[3, k] + tmp

      lω[1, k] = exp(-tmp)
      lω[2, k] = exp(logq[4, k] - logq[3, k] - tmp)
    else
      tmp = log1p(exp(logq[3, k] - logq[4, k]))

      lγ[3] += logq[4, k] + tmp

      lω[1, k] = exp(logq[3, k] - logq[4, k] - tmp)
      lω[2, k] = exp(-tmp)
    end

    if (lγ[1] >= lγ[2]) && (lγ[1] >= lγ[3])
      tmp = lγ[1] + log1p(exp(lγ[2] - lγ[1]) + exp(lγ[3] - lγ[1]))
    elseif (lγ[2] >= lγ[1]) && (lγ[2] >= lγ[3])
      tmp = lγ[2] + log1p(exp(lγ[1] - lγ[2]) + exp(lγ[3] - lγ[2]))
    else
      tmp = lγ[3] + log1p(exp(lγ[1] - lγ[3]) + exp(lγ[2] - lγ[3]))
    end

    u = rand()

    if u <= exp(lγ[1] - tmp)
      for g in 1:k
        C[g, b] = 0x01
      end
      logpr += logγ[1] + logω[1]
      logpo += lγ[1] - tmp
    elseif u <= exp(lγ[1] - tmp) + exp(lγ[2] - tmp)
      for g in 1:k
        C[g, b] = 0x02
      end
      logpr += logγ[2] + logω[2]
      logpo += lγ[2] - tmp
    else
      for g in 1:k
        if rand() <= exp(lω[1, l])
          C[g, b] = 0x03
          logpr += logω[3]
          logpo += lω[1, l]
        else
          C[g, b] = 0x04
          logpr += logω[4]
          logpo += lω[2, l]
        end
      end
      logpr += logγ[3] + logω[3]
      logpo += lγ[3] - tmp
    end
  end

  (logpr, logpo)
end
