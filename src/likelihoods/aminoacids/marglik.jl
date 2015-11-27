# This file is part of Kpax3. License is MIT.

function logmarglik(x::Array{Float64, 1},
                    n::Real,
                    α::Array{Float64, 1},
                    β::Array{Float64, 1})
  lbeta(α + x, β + n - x) - lbeta(α, β)
end

function logmarglik(x::Real,
                    n::Real,
                    α::Float64,
                    β::Float64)
  lbeta(α + x, β + n - x) - lbeta(α, β)
end

function logcondmarglik(x::Array{UInt8, 1},
                        n::Real,
                        n1s::Array{Float64, 1},
                        α::Array{Float64, 1},
                        β::Array{Float64, 1})
  Float64(x) .* log(α + n1s) +
  Float64(0x01 - x) .* log(β + n - n1s) -
  log(α + β + n)
end

function logcondmarglik(x::UInt8,
                        n::Real,
                        n1s::Float64,
                        α::Float64,
                        β::Float64)
  Float64(x) * log(α + n1s) +
  Float64(0x01 - x) * log(β + n - n1s) -
  log(α + β + n)
end

function marglik(x::Array{Float64, 1},
                 n::Real,
                 α::Array{Float64, 1},
                 β::Array{Float64, 1})
  exp(lbeta(α + x, β + n - x) - lbeta(α, β))
end

function marglik(x::Real,
                 n::Real,
                 α::Float64,
                 β::Float64)
  exp(lbeta(α + x, β + n - x) - lbeta(α, β))
end

function condmarglik(x::Array{UInt8, 1},
                     n::Real,
                     n1s::Array{Float64, 1},
                     α::Array{Float64, 1},
                     β::Array{Float64, 1})
  exp(x * log(α + n1s) + (0x01 - x) * log(β + n - n1s) - log(α + β + n))
end

function condmarglik(x::UInt8,
                     n::Real,
                     n1s::Float64,
                     α::Float64,
                     β::Float64)
  exp(x * log(α + n1s) + (0x01 - x) * log(β + n - n1s) - log(α + β + n))
end
