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

function normalizepartition(infile::AbstractString,
                            n::Int)
  d = readcsv(infile, Int)

  if size(d, 2) != 1
    throw(KInputError("Too many columns found in file ", infile, "."))
  end

  if length(d) != n
    throw(KInputError(string("Partition length is not equal to the sample ",
                             "size: ", length(d), " instead of ", n, ".")))
  end

  indexin(d[:, 1], sort(unique(d[:, 1])))
end

function normalizepartition(infile::AbstractString,
                            id::Vector{ASCIIString})
  d = readcsv(infile, ASCIIString)

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
        throw(KInputError(string("Missing ids in file ", infile, ".")))
      end
    end

    indexin(partition, sort(unique(partition)))[idx]
  else
    throw(KInputError("Too many columns found in file ", infile, "."))
  end
end

function initializepartition(D::Matrix{Float64},
                             priorR::PriorRowPartition,
                             priorC::AminoAcidPriorCol;
                             kset::UnitRange{Int}=1:0)
  n = size(D, 1)

  if length(kset) == 0
    kset = 2:max(ceil(Int, sqrt(n)), 100)
  elseif kset[1] > 0
    if kset[end] > n
      kset = (kset[1] == 1) ? 2:n : kset[1]:n
    elseif kset[1] == 1
      kset = 2:kset[end]
    end
  else
    throw(KDomainError("First element of 'kset' is less than one."))
  end

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r, maxclust=kset[end])

  R = ones(Int, n)

  s = AminoAcidState(x.data, R, priorR, priorC, settings)
  slp = s.logpR + s.logpC[1] + s.loglik

  D = zeros(Float64, n, n)
  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    D[i, j] = D[j, i] = d[idx]
    idx += 1
  end

  t1 = copystate(s)
  tlp1 = slp

  t2 = copystate(s)
  tlp2 = slp

  niter = 0
  for k in kset
    updateprior!(priorC, k)

    copy!(R, kmedoids(D, k).assignments)
    updatestate!(t1, x.data, R, priorR, priorC, settings)
    tlp1 = t1.logpR + t1.logpC[1] + t1.loglik

    niter = 0
    while niter < 10
      copy!(R, kmedoids(D, k).assignments)
      updatestate!(t2, x.data, R, priorR, priorC, settings)
      tlp2 = t2.logpR + t2.logpC[1] + t2.loglik

      if tlp2 > tlp1
        copystate!(t1, t2)
        tlp1 = tlp2
      end

      niter += 1
    end

    if tlp1 > slp
      copystate!(s, t1)
      slp = tlp1
    end
  end

end
