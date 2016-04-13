# This file is part of Kpax3. License is MIT.

"""
remove "gaps" and non-positive values from the partition
Example: [1, 1, 0, 1, -2, 4, 0] -> [3, 3, 2, 3, 1, 4, 2]
"""
function initpartition(partition::Vector{Int},
                       n::Int)
  if length(partition) != n
    throw(KInputError(string("Length of argument 'partition' is not equal to ",
                             "the sample size: ", length(partition),
                             " instead of ", n, ".")))
  end

  indexin(partition, sort(unique(partition)))
end

function initpartition(infile::AbstractString,
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

function initpartition(infile::AbstractString,
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
