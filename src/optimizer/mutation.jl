# This file is part of Kpax3. License is MIT.

function mutation!(o::KOffspring,
                   mrate::Float64)
  n = length(o.R)

  w = zeros(Int, n)
  found = false
  i = 0

  for a in 1:n
    if rand() <= mrate
      # remove a from its cluster
      o.v[o.R[a]] -= 1

      copy!(w, o.v)

      # candidate empty cluster to be filled
      found = false
      i = 0
      while (!found) && (i < n)
        i += 1
        if w[i] == 0
          w[i] = 1
          found = true
        end
      end

      o.R[a] = StatsBase.sample(WeightVec(w, n))
      o.v[R[a]] += 1
    end
  end

  nothing
end
