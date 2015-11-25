# This file is part of Kpax3. License is MIT.

type AminoAcidPriorCol <: PriorColPartition
  γ::Array{Float64, 1}
  logγ::Array{Float64, 1}

  ω::Array{Float64, 1}
  logω::Array{Float64, 1}

  A::Array{Float64, 2}
  B::Array{Float64, 2}
end

function AminoAcidPriorCol(data::Array{UInt8, 2},
                           k::Int,
                           g::Array{Float64, 1},
                           r::Float64)
  if length(g) != 3
    throw(Kpax3InputError("Argument 'g' does not have length 3."))
  end

  if any(g .< 0)
    throw(Kpax3DomainError("Argument 'g' contains negative values."))
  end

  # probabilities must sum to one
  g = g / (g[1] + g[2] + g[3])

  m, n = size(data)

  # attributes 3 and 4 are both from property 3 -> use g.prob[3] twice
  γ = [g[1], g[2], g[3], g[3]]

  # log(0) is -Inf
  # since we are going to multiply this matrix with a positive scalar (1 / k)
  # no NaN can be produced even if some γ[i] are zero
  logγ = log(γ)

  ω = [1.0, 1.0, (k - 1.0) / k, 1.0 / k]
  logω = [0.0, 0.0, log(k - 1.0) - log(k), -log(k)]

  A = zeros(Float64, m, 4)
  B = zeros(Float64, m, 4)

  n1s = zeros(Float64, m)
  for j = 1:m
    n1s[j] = data[j, 1]
  end
  for i = 2:n, j = 1:m
    n1s[j] += data[j, i]
  end

  # uninformative attributes
  # If n > r, these two parameters sum to (r+1), i.e. the theoretical sample
  # size for the characteristic hyperparameters. This is done so because we
  # don't want them to be overwhelmed by the data.
  # The mean is the same to the one obtained with a Jeffreys prior
  if n > r
    A[:, 1] = (r + 1.0) * (n1s + 0.5) / (n + 1)
    B[:, 1] = (r + 1.0) - A[:, 1]
  else
    A[:, 1] = n1s + 0.5
    B[:, 1] = n - n1s + 0.5
  end

  # informative but not characteristic for any cluster
  A[:, 2] = 1.0
  B[:, 2] = 1.0

  # informative and characteristic... but for another cluster
  A[:, 3] = 1.0
  B[:, 3] = r

  # informative and characteristic for this cluster
  A[:, 4] = r
  B[:, 4] = 1.0

  AminoAcidPriorCol(γ, logγ, ω, logω, A, B)
end
