# This file is part of Kpax3. License is MIT.

type Cluster
  v::Int
  unit::Array{Int, 1}
  n1s::Array{Float64, 1}
  ll::Float64
end

#=
function Cluster(g::Int,
                 partition::Array{Int, 1},
                 data::Array{UInt8, 2},
                 maxsize::Int)
  m, n = size(data)

  v = 0
  unit = zeros(Int, min(maxsize, n))
  n1s = zeros(Float64, m)

  l = 0

  for i in 1:n
    if partition[i] == g
      v += 1
      l += 1

      if l > maxsize
        maxsize *= 2

        tmp = zeros(Int, maxsize)
        tmp[1:(l - 1)] = unit

        unit = tmp
      end

      unit[l] = i

      for j in 1:m
        n1s[j] += data[j, i]
      end
    end
  end

  Cluster(v, unit, n1s)
end
=#
