# This file is part of Kpax3. License is MIT.

function rcolpartition!(C::Array{UInt8, 2},
                        priorC::AminoAcidPriorCol,
                        cluster::Array{Cluster, 1},
                        emptycluster::BitArray{1})
  cl = find(!emptycluster)
  k = length(cl)

  m = size(C, 1)

  logq = zeros(Float64, 4, k)

  γ = zeros(Float64, 3)
  logγ = zeros(Float64, 3)

  ω = zeros(Float64, 2, k)

  g = zero(Int)

  for j in 1:m
    logγ[1] = priorC.logγ[1]
    logγ[2] = priorC.logγ[2]
    logγ[3] = priorC.logγ[3]

    for i in 1:k
      g = cl[i]

      logq[1, i] = logmarglik(cluster[g].n1s[j], cluster[g].v, priorC.A[1, j],
                              priorC.B[1, j])
      logq[2, i] = logmarglik(cluster[g].n1s[j], cluster[g].v, priorC.A[2, j],
                              priorC.B[2, j])
      logq[3, i] = priorC.logω[3] + logmarglik(cluster[g].n1s[j], cluster[g].v,
                                               priorC.A[3, j], priorC.B[3, j])
      logq[4, i] = priorC.logω[4] + logmarglik(cluster[g].n1s[j], cluster[g].v,
                                               priorC.A[4, j], priorC.B[4, j])

      tmp = log(exp(logq[3, i]) + exp(logq[4, i]))

      logγ[1] += logq[1, i]
      logγ[2] += logq[2, i]
      logγ[3] += tmp

      ω[1, i] = exp(logq[3, i] - tmp)
      ω[2, i] = exp(logq[4, i] - tmp)
    end

    tmp = log(exp(logγ[1]) + exp(logγ[2]) + exp(logγ[3]))

    γ[1] = exp(logγ[1] - tmp)
    γ[2] = exp(logγ[2] - tmp)
    γ[3] = exp(logγ[3] - tmp)

    u = sample(StatsBase.WeightVec(γ))

    if u == 1
      C[j, cl] = 0x01
    elseif u == 2
      C[j, cl] = 0x02
    else
      for i in 1:k
        C[j, cl[i]] = UInt8(2 + sample(WeightVec(vec(ω[:, i]))))
      end
    end
  end

  nothing
end
