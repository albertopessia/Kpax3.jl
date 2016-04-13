# This file is part of Kpax3. License is MIT.

function loss_binder(R::Vector{Int},
                     P::Vector{Float64})
  n = length(R)

  loss = 0.0
  idx = 0
  for j in 1:(n - 1), i in (j + 1):n
    idx += 1
    loss += (R[i] == R[j]) ? 1 - P[idx] : P[idx]
  end

  loss
end
