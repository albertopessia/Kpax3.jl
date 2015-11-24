# This file is part of Kpax3. License is MIT.

function logmarglik(x::Array{Float64, 1},
                    n::Float64,
                    α::Array{Float64, 1},
                    β::Array{Float64, 1})
  lbeta(α + x, β + n - x) - lbeta(α, β)
end

function logmarglik(x::Float64,
                    n::Float64,
                    α::Float64,
                    β::Float64)
  lbeta(α + x, β + n - x) - lbeta(α, β)
end

function marglik(x::Array{Float64, 1},
                 n::Float64,
                 α::Array{Float64, 1},
                 β::Array{Float64, 1})
  exp(lbeta(α + x, β + n - x) - lbeta(α, β))
end

function marglik(x::Float64,
                 n::Float64,
                 α::Float64,
                 β::Float64)
  exp(lbeta(α + x, β + n - x) - lbeta(α, β))
end
