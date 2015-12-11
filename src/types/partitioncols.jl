# This file is part of Kpax3. License is MIT.

abstract PriorColPartition

type AminoAcidPriorCol <: PriorColPartition
  γ::Vector{Float64}
  logγ::Vector{Float64}

  ω::Vector{Float64}
  logω::Vector{Float64}

  A::Matrix{Float64}
  B::Matrix{Float64}
end

function AminoAcidPriorCol(data::Matrix{UInt8},
                           k::Int,
                           g::Vector{Float64},
                           r::Float64)
  if length(g) != 3
    throw(KInputError("Argument 'g' does not have length 3."))
  end

  if g[1] < 0
    throw(KDomainError("Argument 'g[1]' is negative."))
  end

  if g[2] < 0
    throw(KDomainError("Argument 'g[2]' is negative."))
  end

  if g[3] < 0
    throw(KDomainError("Argument 'g[3]' is negative."))
  end

  m, n = size(data)

  # probabilities must sum to one
  tot = g[1] + g[2] + g[3]
  g[1] /= tot
  g[2] /= tot
  g[3] /= tot

  # attributes 3 and 4 are both from property 3 -> use g[3] twice
  γ = [g[1]; g[2]; g[3]; g[3]]

  # log(0) is -Inf
  # since we are going to multiply this matrix with a positive scalar (1 / k)
  # no NaN can be produced even if some γ[i] are zero
  logγ = [log(γ[1]); log(γ[2]); log(γ[3]); log(γ[3])]

  ω = [1.0; 1.0; (k - 1.0) / k; 1.0 / k]
  logω = [0.0; 0.0; log(k - 1.0) - log(k); -log(k)]

  n1s = zeros(Float64, m)
  for a in 1:n, b in 1:m
    n1s[b] += Float64(data[b, a])
  end

  A = zeros(Float64, 4, m)
  B = zeros(Float64, 4, m)

  if n > r
    for b in 1:m
      # uninformative attributes
      # If n > r, these two parameters sum to (r+1), i.e. the theoretical sample
      # size for the characteristic hyperparameters. This is done so because we
      # don't want them to be overwhelmed by the data.
      # The mean is the same to the one obtained with a Jeffreys prior
      A[1, b] = (r + 1.0) * (n1s[b] + 0.5) / (n + 1.0)
      B[1, b] = (r + 1.0) - A[1, b]

      # informative but not characteristic for any cluster
      A[2, b] = 1.0
      B[2, b] = 1.0

      # informative and characteristic... but for another cluster
      A[3, b] = 1.0
      B[3, b] = r

      # informative and characteristic for this cluster
      A[4, b] = r
      B[4, b] = 1.0
    end
  else
    for b in 1:m
      # uninformative attributes
      A[1, b] = n1s[b] + 0.5
      B[1, b] = n - n1s[b] + 0.5

      # informative but not characteristic for any cluster
      A[2, b] = 1.0
      B[2, b] = 1.0

      # informative and characteristic... but for another cluster
      A[3, b] = 1.0
      B[3, b] = r

      # informative and characteristic for this cluster
      A[4, b] = r
      B[4, b] = 1.0
    end
  end

  AminoAcidPriorCol(γ, logγ, ω, logω, A, B)
end
