# This file is part of Kpax3. License is MIT.

"""
logp(S | R, X)
"""
function logcondpostS(S::Array{UInt8, 1},
                      priorC::AminoAcidPriorCol,
                      cluster::Array{Cluster, 1},
                      emptycluster::BitArray{1})
  cl = find(!emptycluster)
  k = length(cl)

  m = size(C, 1)

  logp = zero(Float64)

  logγ = zeros(Float64, 3)

  g = zero(Int)

  for j in 1:m
    logγ[1] = priorC.logγ[1]
    logγ[2] = priorC.logγ[2]
    logγ[3] = priorC.logγ[3]

    for i in 1:k
      g = cl[i]

      logγ[1] += logmarglik(cluster[g].n1s[j], cluster[g].v, priorC.A[j, 1],
                            priorC.B[j, 1])
      logγ[2] += logmarglik(cluster[g].n1s[j], cluster[g].v, priorC.A[j, 2],
                            priorC.B[j, 2])
      logγ[3] += log(exp(priorC.logω[3] + logmarglik(cluster[g].n1s[j],
                                                     cluster[g].v,
                                                     priorC.A[j, 3],
                                                     priorC.B[j, 3])) +
                     exp(priorC.logω[4] + logmarglik(cluster[g].n1s[j],
                                                     cluster[g].v,
                                                     priorC.A[j, 4],
                                                     priorC.B[j, 4])))
    end

    logp += (logγ[S[j]] - log(exp(logγ[1]) + exp(logγ[2]) + exp(logγ[3])))
  end

  logp
end

"""
Logp(C | R, X)
"""
function logcondpostC(C::Array{UInt8, 2},
                      priorC::AminoAcidPriorCol,
                      cluster::Array{Cluster, 1},
                      emptycluster::BitArray{1})
  cl = find(!emptycluster)
  k = length(cl)

  m = size(C, 1)

  logp = zero(Float64)

  logq = zeros(Float64, 4, k)
  logγ = zeros(Float64, 3)

  tmp = zero(Float64)
  g = zero(Int)

  for j in 1:m
    logγ[1] = priorC.logγ[1]
    logγ[2] = priorC.logγ[2]
    logγ[3] = priorC.logγ[3]

    for i in 1:k
      g = cl[i]

      logq[1, i] = logmarglik(cluster[g].n1s[j], cluster[g].v, priorC.A[j, 1],
                              priorC.B[j, 1])
      logq[2, i] = logmarglik(cluster[g].n1s[j], cluster[g].v, priorC.A[j, 2],
                              priorC.B[j, 2])
      logq[3, i] = priorC.logω[3] + logmarglik(cluster[g].n1s[j], cluster[g].v,
                                               priorC.A[j, 3], priorC.B[j, 3])
      logq[4, i] = priorC.logω[4] + logmarglik(cluster[g].n1s[j], cluster[g].v,
                                               priorC.A[j, 4], priorC.B[j, 4])

      tmp = log(exp(logq[3, i]) + exp(logq[4, i]))

      logγ[1] += logq[1, i]
      logγ[2] += logq[2, i]
      logγ[3] += tmp

      if C[j, g] == 3
        logp += logq[3, i] - tmp
      elseif C[j, g] == 4
        logp += logq[4, i] - tmp
      end
    end

    if C[j, cl[1]] == 1
      logp += logγ[1] - log(exp(logγ[1]) + exp(logγ[2]) + exp(logγ[3]))
    elseif C[j, cl[1]] == 2
      logp += logγ[2] - log(exp(logγ[1]) + exp(logγ[2]) + exp(logγ[3]))
    else
      logp += logγ[3] - log(exp(logγ[1]) + exp(logγ[2]) + exp(logγ[3]))
    end
  end

  logp
end
