# This file is part of Kpax3. License is MIT.

function rcolpartition!(x::AminoAcidPriorCol,
                        C::Array{UInt8, 2},
                        v::Array{Float64, 1},
                        n1s::Array{Float64, 2},
                        emptyclusters::BitArray{1})
  cl = find(!emptyclusters)

  m = size(C, 1)
  k = length(cl)

  w = zeros(Float64, m, k, 3)
  p = zeros(Float64, m, k, 2)

  # for the moment, use logarithms for a numerically stable algorithm
  # at the end we will proceed to exp them
  γ = zeros(Float64, m, 3)

  tmp = zeros(Float64, 4)
  s = zero(Float64)
  u = zero(Int)

  for i in 1:k, j in 1:m
    g = cl[i]

    tmp[1] = marglik(n1s[j, g], v[g], x.A[j, 1], x.B[j, 1])
    tmp[2] = marglik(n1s[j, g], v[g], x.A[j, 2], x.B[j, 2])
    tmp[3] = x.ω[3] * marglik(n1s[j, g], v[g], x.A[j, 3], x.B[j, 3])
    tmp[4] = x.ω[4] * marglik(n1s[j, g], v[g], x.A[j, 4], x.B[j, 4])

    w[j, i, 1] = tmp[1]
    w[j, i, 2] = tmp[2]
    w[j, i, 3] = tmp[3] + tmp[4]

    p[j, i, 1] = tmp[3] / w[j, i, 3]
    p[j, i, 2] = tmp[4] / w[j, i, 3]

    γ[j, 1] += log(w[j, i, 1])
    γ[j, 2] += log(w[j, i, 2])
    γ[j, 3] += log(w[j, i, 3])
  end

  for j in 1:m
    tmp[1] = exp(log(x.γ[1]) + γ[j, 1])
    tmp[2] = exp(log(x.γ[2]) + γ[j, 2])
    tmp[3] = exp(log(x.γ[3]) + γ[j, 3])

    s = tmp[1] + tmp[2] + tmp[3]

    γ[j, 1] = tmp[1] / s
    γ[j, 2] = tmp[2] / s
    γ[j, 3] = tmp[3] / s
  end

  for j in 1:m
    u = sample(WeightVec(vec(γ[j, :])))

    if u == 1
      C[j, cl] = 0x01
    elseif u == 2
      C[j, cl] = 0x02
    else
      for i in 1:k
        C[j, cl[i]] = UInt8(2 + sample(WeightVec(vec(p[j, i, :]))))
      end
    end
  end

  nothing
end
