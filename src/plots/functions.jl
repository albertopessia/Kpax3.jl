# This file is part of Kpax3. License is MIT.

function reorderunits(R::Vector{Int},
                      P::Matrix{Float64},
                      clusterorder::Vector{Int})
  n = length(R)
  k = maximum(R)

  M = zeros(Float64, n)
  v = zeros(Float64, n)
  for j in 1:(n - 1), i in (j + 1):n
    if R[i] == R[j]
      M[i] += P[i, j]
      v[i] += 1

      M[j] += P[i, j]
      v[j] += 1
    end
  end

  M ./= v

  neworder = zeros(Int, n)
  midpoint = zeros(Float64, k)
  seppoint = zeros(Float64, k + 1)

  h = 1
  u = 1
  for g in 1:k
    idx = find(R .== clusterorder[g])
    u = length(idx)
    ord = sortperm(M[idx], rev=true)

    copy!(neworder, h, idx[ord], 1, u)
    midpoint[g] = (2 * h + u - 1) / 2
    seppoint[g] = h - 0.5

    h += u
  end

  seppoint[k + 1] = n + 0.5

  (neworder, midpoint, seppoint)
end

function computefontsize(width::Real,
                         height::Real)
  # coefficients estimated with a linear model. They can be improved
  Measures.Length(:mm, 1.8933 + 0.0039 * width + 0.0141 * height)
end
