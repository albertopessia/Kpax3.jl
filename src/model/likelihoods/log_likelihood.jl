# This file is part of Kpax3. License is MIT.

function loglikelihood(C::Matrix{UInt8},
                       cl::Vector{Int},
                       k::Int,
                       v::Vector{Int},
                       n1s::Matrix{Float64},
                       priorC::AminoAcidPriorCol)
  m = size(C, 2)
  loglik = 0.0

  # If array A has dimension (d_{1}, ..., d_{l}, ..., d_{L}), to access
  # element A[i_{1}, ..., i_{l}, ..., i_{L}] it is possible to use the
  # following linear index
  # lidx = i_{1} + d_{1} * (i_{2} - 1) + ... +
  #      + (d_{1} * ... * d_{l-1}) * (i_{l} - 1) + ... +
  #      + (d_{1} * ... * d_{L-1}) * (i_{L} - 1)
  #
  # A[i_{1}, ..., i_{l}, ..., i_{L}] == A[lidx]
  lidx = 0
  for b in 1:m
    for l in 1:k
      lidx = C[cl[l], b] + 4 * (b - 1)
      loglik += logmarglik(n1s[cl[l], b], v[cl[l]], priorC.A[lidx],
                           priorC.B[lidx])
    end
  end

  loglik
end
