# This file is part of Kpax3. License is MIT.

"""
remove "gaps" and non-positive values from the partition
Example: [1, 1, 0, 1, -2, 4, 0] -> [3, 3, 2, 3, 1, 4, 2]
"""
function normalizepartition(partition::Vector{Int},
                            n::Int)
  if length(partition) != n
    throw(KInputError(string("Length of argument 'partition' is not equal to ",
                             "the sample size: ", length(partition),
                             " instead of ", n, ".")))
  end

  indexin(partition, sort(unique(partition)))
end

function normalizepartition(ifile::AbstractString,
                            n::Int)
  d = readcsv(ifile, Int)

  if size(d, 2) != 1
    throw(KInputError("Too many columns found in file ", ifile, "."))
  end

  if length(d) != n
    throw(KInputError(string("Partition length is not equal to the sample ",
                             "size: ", length(d), " instead of ", n, ".")))
  end

  indexin(d[:, 1], sort(unique(d[:, 1])))
end

function normalizepartition(ifile::AbstractString,
                            id::Vector{ASCIIString})
  d = readcsv(ifile, ASCIIString)

  if size(d, 1) != length(id)
    throw(KInputError(string("Partition length is not equal to the sample ",
                             "size: ", size(d, 1), " instead of ", length(id),
                             ".")))
  end

  if size(d, 2) == 1
    partition = [parse(Int, x) for x in d[:, 1]]
    indexin(partition, sort(unique(partition)))
  elseif size(d, 2) == 2
    idx = indexin([x::ASCIIString for x in d[:, 1]], id)
    partition = [parse(Int, x) for x in d[:, 2]]

    for i in 1:length(idx)
      if idx[i] == 0
        throw(KInputError(string("Missing ids in file ", ifile, ".")))
      end
    end

    indexin(partition, sort(unique(partition)))[idx]
  else
    throw(KInputError("Too many columns found in file ", ifile, "."))
  end
end

function encodepartition(partition::Vector{Int})
  n = length(partition)

  R = zeros(Int, n)
  lidx = zeros(Int, n)

  g = 0
  for a in 1:n
    g = partition[a]

    if lidx[g] == 0
      lidx[g] = a
    end

    # linear index of unit a in the n-by-n adjacency matrix
    R[a] = a + n * (lidx[g] - 1)
  end

  R
end

function encodepartition!(R::Vector{Int})
  n = length(R)

  lidx = zeros(Int, n)

  g = 0
  for a in 1:n
    g = R[a]

    if lidx[g] == 0
      lidx[g] = a
    end

    # linear index of unit a in the n-by-n adjacency matrix
    R[a] = a + n * (lidx[g] - 1)
  end

  nothing
end

function modifypartition!(R::Vector{Int},
                          k::Int)
  n = length(R)
  q = k

  if (k > 0) && (k < n + 1)
    q = rand(max(1, k - 10):min(n, k + 10))

    if q < k
      modifymerge!(R, k, q)
    elseif q > k
      modifysplit!(R, k, q)
    else
      modifyscramble!(R, k)
    end
  end

  q
end

function modifymerge!(R::Vector{Int},
                      k::Int,
                      q::Int)
  n = length(R)

  cset = zeros(Int, k)
  empty = trues(n)
  l = 1
  for a in 1:n
    if empty[R[a]]
      empty[R[a]] = false
      cset[l] = R[a]
      l += 1
    end
  end

  c = zeros(Int, 2)

  while k != q
    StatsBase.sample!(cset, c, replace=false, ordered=false)

    for a in 1:n
      if R[a] == c[2]
        R[a] = c[1]
      end
    end

    l = 1
    while cset[l] != c[2]
      l += 1
    end

    deleteat!(cset, l)
    k -= 1
  end

  nothing
end

function modifysplit!(R::Vector{Int},
                      k::Int,
                      q::Int)
  n = length(R)

  t = zeros(Int, n)
  for a in 1:n
    t[R[a]] += 1
  end

  g = 0

  while k != q
    k += 1

    w = StatsBase.WeightVec(Float64[t[a] > 0 ? t[a] - 1 : 0 for a in 1:n])
    g = StatsBase.sample(w)

    for a in 1:n
      if (R[a] == g) && ((t[k] == 0) || (rand() <= 0.25))
        R[a] = k
        t[g] -= 1
        t[k] += 1
      end
    end
  end

  nothing
end

function modifyscramble!(R::Vector{Int},
                         k::Int)
  n = length(R)

  t = zeros(Int, n)
  for a in 1:n
    t[R[a]] += 1
  end

  v = zeros(Int, n)
  moved = false

  l = 1
  g = 1
  h = 1
  while g < k
    if t[l] > 0
      t[l] = 0
      v[l] = 0

      for a in (l + 1):n
        v[a] = t[a] > 0 ? t[a] - 1 : 0
      end

      w = StatsBase.WeightVec(v)
      h = StatsBase.sample(w)

      keepgoing = true
      moved = false
      a = 1
      while keepgoing
        if R[a] == h
          if moved
            if rand() <= 0.05
              R[a] = g
              t[h] -= 1

              if t[h] == 1
                keepgoing = false
              end
            end
          else
            moved = true
            R[a] = g
            t[h] -= 1

            if t[h] == 1
              keepgoing = false
            end
          end
        end

        a += 1
        if a > n
          keepgoing = false
        end
      end

      g += 1
    end

    l += 1
  end

  nothing
end
