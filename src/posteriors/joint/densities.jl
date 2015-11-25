# This file is part of Kpax3. License is MIT.

"""
Logposterior minus the log(normalizing constant)
"""
function logposterior(R::Array{Int, 1},
                      C::Array{UInt8, 2},
                      priorR::PriorRowPartition,
                      priorC::PriorColPartition,
                      v::Array{Float64, 1},
                      n1s::Array{Int, 2},
                      emptyclusters::BitArray{1})
  cl = find(!emptyclusters)

  m = size(C, 1)
  k = length(cl)

  logp = logdPriorRow(priorR, R)

  g = cl[1]
  s = zero(UInt8)

  for j in 1:m
    s = C[j, g]
    logp += (priorC.logγ[s] + priorC.logω[s] +
             logmarglik(n1s[j, g], v[g], priorC.A[j, s], priorC.B[j, s]))
  end

  for i in 2:k
    g = cl[i]

    for j in 1:m
      s = C[j, g]
      logp += (priorC.logω[s] +
               logmarglik(n1s[j, g], v[g], priorC.A[j, s], priorC.B[j, s]))
    end
  end

  logp
end
