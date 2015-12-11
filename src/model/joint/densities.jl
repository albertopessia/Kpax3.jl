# This file is part of Kpax3. License is MIT.

"""
Logposterior minus the log(normalizing constant)
"""
function logposterior(R::Vector{Int},
                      C::Matrix{UInt8},
                      priorR::PriorRowPartition,
                      priorC::PriorColPartition,
                      k::Int)
#=
  logp = logdPriorRow(priorR, R)
  s = zero(UInt8)

  for b in 1:size(C, 2)
    s = C[1, b]
    logp += (priorC.logγ[s] + priorC.logω[s] +
             logmarglik(cluster[1].n1s[b], cluster[1].v, priorC.A[s, b],
                        priorC.B[s, b]))

    for g in 2:k
      s = C[g, b]
      logp += (priorC.logω[s] + logmarglik(cluster[g].n1s[b], cluster[g].v,
                                           priorC.A[s, b], priorC.B[s, b]))
    end
  end

  logp
=#
end
