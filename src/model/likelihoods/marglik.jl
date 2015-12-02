# This file is part of Kpax3. License is MIT.

function logmarglik(y::Array{Float64, 1},
                    n::Real,
                    α::Array{Float64, 1},
                    β::Array{Float64, 1})
  lbeta(α + y, β + n - y) - lbeta(α, β)
end

function logmarglik(y::Real,
                    n::Real,
                    α::Float64,
                    β::Float64)
  lbeta(α + y, β + n - y) - lbeta(α, β)
end

function logcondmarglik(x::Array{UInt8, 1},
                        y::Array{Float64, 1},
                        n::Real,
                        α::Array{Float64, 1},
                        β::Array{Float64, 1})
  [(x[i] * log(α[i] + y[i]) + (0x01 - x[i]) * log(β[i] + n - y[i]) -
    log(α[i] + β[i] + n))::Float64 for i in 1:length(x)]
end

function logcondmarglik(x::UInt8,
                        y::Float64,
                        n::Real,
                        α::Float64,
                        β::Float64)
  (x * log(α + y) + (0x01 - x) * log(β + n - y) - log(α + β + n))::Float64
end

function marglik(y::Array{Float64, 1},
                 n::Real,
                 α::Array{Float64, 1},
                 β::Array{Float64, 1})
  exp(lbeta(α + y, β + n - y) - lbeta(α, β))
end

function marglik(y::Real,
                 n::Real,
                 α::Float64,
                 β::Float64)
  exp(lbeta(α + y, β + n - y) - lbeta(α, β))
end

function condmarglik(x::Array{UInt8, 1},
                     y::Array{Float64, 1},
                     n::Real,
                     α::Array{Float64, 1},
                     β::Array{Float64, 1})
  [exp(x[i] * log(α[i] + y[i]) + (0x01 - x[i]) * log(β[i] + n - y[i]) -
       log(α[i] + β[i] + n))::Float64 for i in 1:length(x)]
end

function condmarglik(x::UInt8,
                     y::Float64,
                     n::Real,
                     α::Float64,
                     β::Float64)
  exp(x * log(α + y) + (0x01 - x) * log(β + n - y) - log(α + β + n))::Float64
end
