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

function initializestates(x::AminoAcidData,
                          d::Vector{Float64},
                          N::Int,
                          settings::KSettings;
                          kset::UnitRange{Int}=1:0)
  m, n = size(x.data)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, 1, settings.γ, settings.r)

  R = ones(Int, n)

  state = AminoAcidState(x.data, R, priorR, priorC, settings)
  logpp = state.logpR + state.logpC[1] + state.loglik

  if length(kset) == 0
    kset = 2:max(ceil(Int, sqrt(n)), 100)
  end

  D = zeros(Float64, n, n)
  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    D[i, j] = D[j, i] = d[idx]
    idx += 1
  end

  tmpstate = copystate(state)
  tmplogpp = logpp

  for k in kset
    updateprior!(priorC, k)

    copy!(R, normalizepartition(kmedoids(D, k; maxiter=1000).assignments, n))

    tmpstate = AminoAcidState(x.data, R, priorR, priorC, settings)
    tmplogpp = tmpstate.logpR + tmpstate.logpC[1] + tmpstate.loglik
  end

    statelist = [copystate(state) for i in 1:N]


end
