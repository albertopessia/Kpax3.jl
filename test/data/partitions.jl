# This file is part of K-Pax3. License is MIT.

"""
# Partitions of n

All the partitions of the set N = {1, ..., n} and their properties.

## Fields

* `n` total number of units (integer to partition)
* `B` total number of partitions of N (Bell number of `n`)
* `C` total number of configuration classes (see details)
* `k` total number of blocks within a particular configuration class
* `partition` explicit partitions of the `n` units
* `blocksize` block sizes corresponding to configuration classes (see details)
* `blockcount` numbers of terms of various sizes (see details)
* `index` indices of configuration classes within variable partition

## Details

Let `k` represent the total number of blocks (clusters).

Configuration classes are the number of distinct ways to partition the integer
`n` into blocks of (n1, ..., nk) sizes, where

* ng = #{i : i ∈ g},
* n1 ≥ ... ng ≥ ... ≥ nk > 0
* n1 + ... + nk = n.

Partitions of `n` can also be represented by the numbers of terms of various
sizes mj = #{g : ng = j} (j = 1, ..., n), where

* m1 + ... mj + ... + mn = k
* 1 * m1 + ... + j * mj + ... + n * mn = n
"""
immutable TestPartition
  n::Int
  B::Int
  C::Int
  k::Array{Int, 1}
  index::Array{Int, 1}
  blocksize::Array{Int, 2}
  blockcount::Array{Int, 2}
  partition::Array{Int, 2}
end

"""Partitions of 1"""
function TestP1()
  n = 1
  B = 1
  C = 1

  k = ones(Int, C)
  index = ones(Int, C)

  blocksize  = ones(Int, n, C)
  blockcount = ones(Int, n, C)

  partition  = ones(Int, n, B)

  TestPartition(n, B, C, k, index, blocksize, blockcount, partition)
end

"""Partitions of 2"""
function TestP2()
  n = 2
  B = 2
  C = 2

  k = zeros(Int, C)
  index = zeros(Int, C)

  blocksize  = zeros(Int, n, C)
  blockcount = zeros(Int, n, C)

  partition  = zeros(Int, n, B)

  i = 0
  j = 0

  # k = 1

  # **
  j += 1
  index[j] = i + 1
  k[j] = 1
  blocksize[:,  j] = [2, 0]
  blockcount[:, j] = [0, 1]

  partition[:, i += 1] = [1, 1]

  # k = 2

  # *|*
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [1, 1]
  blockcount[:, j] = [2, 0]

  partition[:, i += 1] = [1, 2]

  TestPartition(n, B, C, k, index, blocksize, blockcount, partition)
end

"""Partitions of 3"""
function TestP3()
  n = 3
  B = 5
  C = 3

  k = zeros(Int, C)
  index = zeros(Int, C)

  blocksize  = zeros(Int, n, C)
  blockcount = zeros(Int, n, C)

  partition  = zeros(Int, n, B)

  i = 0
  j = 0

  # k = 1

  # ***
  j += 1
  index[j] = i + 1
  k[j] = 1
  blocksize[:,  j] = [3, 0, 0]
  blockcount[:, j] = [0, 0, 1]

  partition[:, i += 1] = [1, 1, 1]

  # k = 2

  # **|*
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [2, 1, 0]
  blockcount[:, j] = [1, 1, 0]

  partition[:, i += 1] = [1, 1, 2]
  partition[:, i += 1] = [1, 2, 1]
  partition[:, i += 1] = [2, 1, 1]

  # k = 3

  # *|*|*
  j += 1
  index[j] = i + 1
  k[j] = 3
  blocksize[:,  j] = [1, 1, 1]
  blockcount[:, j] = [3, 0, 0]

  partition[:, i += 1] = [1, 2, 3]

  TestPartition(n, B, C, k, index, blocksize, blockcount, partition)
end

"""Partitions of 4"""
function TestP4()
  n = 4
  B = 15
  C = 5

  k = zeros(Int, C)
  index = zeros(Int, C)

  blocksize  = zeros(Int, n, C)
  blockcount = zeros(Int, n, C)

  partition  = zeros(Int, n, B)

  i = 0
  j = 0

  # k = 1

  # ****
  j += 1
  index[j] = i + 1
  k[j] = 1
  blocksize[:,  j] = [4, 0, 0, 0]
  blockcount[:, j] = [0, 0, 0, 1]

  partition[:, i += 1] = [1, 1, 1, 1]

  # k = 2

  # ***|*
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [3, 1, 0, 0]
  blockcount[:, j] = [1, 0, 1, 0]

  partition[:, i += 1] = [1, 1, 1, 2]
  partition[:, i += 1] = [1, 1, 2, 1]
  partition[:, i += 1] = [1, 2, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1]

  # **|**
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [2, 2, 0, 0]
  blockcount[:, j] = [0, 2, 0, 0]

  partition[:, i += 1] = [1, 1, 2, 2]
  partition[:, i += 1] = [1, 2, 1, 2]
  partition[:, i += 1] = [1, 2, 2, 1]

  # k = 3

  # **|*|*
  j += 1
  index[j] = i + 1
  k[j] = 3
  blocksize[:,  j] = [2, 1, 1, 0]
  blockcount[:, j] = [2, 1, 0, 0]

  partition[:, i += 1] = [1, 1, 2, 3]
  partition[:, i += 1] = [1, 2, 1, 3]
  partition[:, i += 1] = [1, 2, 3, 1]
  partition[:, i += 1] = [2, 1, 1, 3]
  partition[:, i += 1] = [2, 1, 3, 1]
  partition[:, i += 1] = [2, 3, 1, 1]

  # k = 4

  # *|*|*|*
  j += 1
  index[j] = i + 1
  k[j] = 4
  blocksize[:,  j] = [1, 1, 1, 1]
  blockcount[:, j] = [4, 0, 0, 0]

  partition[:, i += 1] = [1, 2, 3, 4]

  TestPartition(n, B, C, k, index, blocksize, blockcount, partition)
end

"""Partitions of 5"""
function TestP5()
  n = 5
  B = 52
  C = 7

  k = zeros(Int, C)
  index = zeros(Int, C)

  blocksize  = zeros(Int, n, C)
  blockcount = zeros(Int, n, C)

  partition  = zeros(Int, n, B)

  i = 0
  j = 0

  # k = 1

  # *****
  j += 1
  index[j] = i + 1
  k[j] = 1
  blocksize[:,  j] = [5, 0, 0, 0, 0]
  blockcount[:, j] = [0, 0, 0, 0, 1]

  partition[:, i += 1] = [1, 1, 1, 1, 1]

  # k = 2

  # ****|*
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [4, 1, 0, 0, 0]
  blockcount[:, j] = [1, 0, 0, 1, 0]

  partition[:, i += 1] = [1, 1, 1, 1, 2]
  partition[:, i += 1] = [1, 1, 1, 2, 1]
  partition[:, i += 1] = [1, 1, 2, 1, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1, 1]

  # ***|**
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [3, 2, 0, 0, 0]
  blockcount[:, j] = [0, 1, 1, 0, 0]

  partition[:, i += 1] = [1, 1, 1, 2, 2]
  partition[:, i += 1] = [1, 1, 2, 1, 2]
  partition[:, i += 1] = [1, 2, 1, 1, 2]
  partition[:, i += 1] = [2, 1, 1, 1, 2]
  partition[:, i += 1] = [1, 1, 2, 2, 1]
  partition[:, i += 1] = [1, 2, 1, 2, 1]
  partition[:, i += 1] = [2, 1, 1, 2, 1]
  partition[:, i += 1] = [1, 2, 2, 1, 1]
  partition[:, i += 1] = [2, 1, 2, 1, 1]
  partition[:, i += 1] = [2, 2, 1, 1, 1]

  # k = 3

  # **|**|*
  j += 1
  index[j] = i + 1
  k[j] = 3
  blocksize[:,  j] = [2, 2, 1, 0, 0]
  blockcount[:, j] = [1, 2, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 2, 2, 3]
  partition[:, i += 1] = [1, 1, 2, 3, 2]
  partition[:, i += 1] = [1, 1, 3, 2, 2]
  partition[:, i += 1] = [1, 2, 1, 2, 3]
  partition[:, i += 1] = [1, 2, 1, 3, 2]
  partition[:, i += 1] = [1, 3, 1, 2, 2]
  partition[:, i += 1] = [1, 2, 2, 1, 3]
  partition[:, i += 1] = [1, 2, 3, 1, 2]
  partition[:, i += 1] = [1, 3, 2, 1, 2]
  partition[:, i += 1] = [1, 2, 2, 3, 1]
  partition[:, i += 1] = [1, 2, 3, 2, 1]
  partition[:, i += 1] = [1, 3, 2, 2, 1]
  partition[:, i += 1] = [3, 1, 1, 2, 2]
  partition[:, i += 1] = [3, 1, 2, 1, 2]
  partition[:, i += 1] = [3, 1, 2, 2, 1]

  # ***|*|*
  j += 1
  index[j] = i + 1
  k[j] = 3
  blocksize[:,  j] = [3, 1, 1, 0, 0]
  blockcount[:, j] = [2, 0, 1, 0, 0]

  partition[:, i += 1] = [1, 1, 1, 2, 3]
  partition[:, i += 1] = [1, 1, 2, 1, 3]
  partition[:, i += 1] = [1, 1, 2, 3, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 3]
  partition[:, i += 1] = [1, 2, 1, 3, 1]
  partition[:, i += 1] = [1, 2, 3, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1, 3]
  partition[:, i += 1] = [2, 1, 1, 3, 1]
  partition[:, i += 1] = [2, 1, 3, 1, 1]
  partition[:, i += 1] = [2, 3, 1, 1, 1]

  # k = 4

  # **|*|*|*
  j += 1
  index[j] = i + 1
  k[j] = 4
  blocksize[:,  j] = [2, 1, 1, 1, 0]
  blockcount[:, j] = [3, 1, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 2, 3, 4]
  partition[:, i += 1] = [1, 2, 1, 3, 4]
  partition[:, i += 1] = [1, 2, 3, 1, 4]
  partition[:, i += 1] = [1, 2, 3, 4, 1]
  partition[:, i += 1] = [2, 1, 1, 3, 4]
  partition[:, i += 1] = [2, 1, 3, 1, 4]
  partition[:, i += 1] = [2, 1, 3, 4, 1]
  partition[:, i += 1] = [2, 3, 1, 1, 4]
  partition[:, i += 1] = [2, 3, 1, 4, 1]
  partition[:, i += 1] = [2, 3, 4, 1, 1]

  # k = 5

  # *|*|*|*|*
  j += 1
  index[j] = i + 1
  k[j] = 5
  blocksize[:,  j] = [1, 1, 1, 1, 1]
  blockcount[:, j] = [5, 0, 0, 0, 0]

  partition[:, i += 1] = [1, 2, 3, 4, 5]

  TestPartition(n, B, C, k, index, blocksize, blockcount, partition)
end

"""Partitions of 6"""
function TestP6()
  n = 6
  B = 203
  C = 11

  k = zeros(Int, C)
  index = zeros(Int, C)

  blocksize  = zeros(Int, n, C)
  blockcount = zeros(Int, n, C)

  partition  = zeros(Int, n, B)

  i = 0
  j = 0

  # k = 1

  # ******
  j += 1
  index[j] = i + 1
  k[j] = 1
  blocksize[:,  j] = [6, 0, 0, 0, 0, 0]
  blockcount[:, j] = [0, 0, 0, 0, 0, 6]

  partition[:, i += 1] = [1, 1, 1, 1, 1, 1]

  # k = 2

  # *****|*
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [5, 1, 0, 0, 0, 0]
  blockcount[:, j] = [1, 0, 0, 0, 1, 0]

  partition[:, i += 1] = [1, 1, 1, 1, 1, 2]
  partition[:, i += 1] = [1, 1, 1, 1, 2, 1]
  partition[:, i += 1] = [1, 1, 1, 2, 1, 1]
  partition[:, i += 1] = [1, 1, 2, 1, 1, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1, 1, 1]

  # ****|**
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [4, 2, 0, 0, 0, 0]
  blockcount[:, j] = [0, 1, 0, 1, 0, 0]

  partition[:, i += 1] = [1, 1, 1, 1, 2, 2]
  partition[:, i += 1] = [1, 1, 1, 2, 1, 2]
  partition[:, i += 1] = [1, 1, 1, 2, 2, 1]
  partition[:, i += 1] = [1, 1, 2, 1, 1, 2]
  partition[:, i += 1] = [1, 1, 2, 1, 2, 1]
  partition[:, i += 1] = [1, 1, 2, 2, 1, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 1, 2]
  partition[:, i += 1] = [1, 2, 1, 1, 2, 1]
  partition[:, i += 1] = [1, 2, 1, 2, 1, 1]
  partition[:, i += 1] = [1, 2, 2, 1, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1, 1, 2]
  partition[:, i += 1] = [2, 1, 1, 1, 2, 1]
  partition[:, i += 1] = [2, 1, 1, 2, 1, 1]
  partition[:, i += 1] = [2, 1, 2, 1, 1, 1]
  partition[:, i += 1] = [2, 2, 1, 1, 1, 1]

  # ***|***
  j += 1
  index[j] = i + 1
  k[j] = 2
  blocksize[:,  j] = [3, 3, 0, 0, 0, 0]
  blockcount[:, j] = [0, 0, 2, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 1, 2, 2, 2]
  partition[:, i += 1] = [1, 1, 2, 1, 2, 2]
  partition[:, i += 1] = [1, 1, 2, 2, 1, 2]
  partition[:, i += 1] = [1, 1, 2, 2, 2, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 2, 2]
  partition[:, i += 1] = [1, 2, 1, 2, 1, 2]
  partition[:, i += 1] = [1, 2, 1, 2, 2, 1]
  partition[:, i += 1] = [1, 2, 2, 1, 1, 2]
  partition[:, i += 1] = [1, 2, 2, 1, 2, 1]
  partition[:, i += 1] = [1, 2, 2, 2, 1, 1]

  # k = 3

  # ****|*|*
  j += 1
  index[j] = i + 1
  k[j] = 3
  blocksize[:,  j] = [4, 1, 1, 0, 0, 0]
  blockcount[:, j] = [2, 0, 0, 1, 0, 0]

  partition[:, i += 1] = [1, 1, 1, 1, 2, 3]
  partition[:, i += 1] = [1, 1, 1, 2, 1, 3]
  partition[:, i += 1] = [1, 1, 1, 2, 3, 1]
  partition[:, i += 1] = [1, 1, 2, 1, 1, 3]
  partition[:, i += 1] = [1, 1, 2, 1, 3, 1]
  partition[:, i += 1] = [1, 1, 2, 3, 1, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 1, 3]
  partition[:, i += 1] = [1, 2, 1, 1, 3, 1]
  partition[:, i += 1] = [1, 2, 1, 3, 1, 1]
  partition[:, i += 1] = [1, 2, 3, 1, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1, 1, 3]
  partition[:, i += 1] = [2, 1, 1, 1, 3, 1]
  partition[:, i += 1] = [2, 1, 1, 3, 1, 1]
  partition[:, i += 1] = [2, 1, 3, 1, 1, 1]
  partition[:, i += 1] = [2, 3, 1, 1, 1, 1]

  # ***|**|*
  j += 1
  index[j] = i + 1
  k[j] = 3
  blocksize[:,  j] = [3, 2, 1, 0, 0, 0]
  blockcount[:, j] = [1, 1, 1, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 1, 2, 2, 3]
  partition[:, i += 1] = [1, 1, 1, 2, 3, 2]
  partition[:, i += 1] = [1, 1, 1, 3, 2, 2]
  partition[:, i += 1] = [1, 1, 2, 1, 2, 3]
  partition[:, i += 1] = [1, 1, 2, 1, 3, 2]
  partition[:, i += 1] = [1, 1, 3, 1, 2, 2]
  partition[:, i += 1] = [1, 1, 2, 2, 1, 3]
  partition[:, i += 1] = [1, 1, 2, 3, 1, 2]
  partition[:, i += 1] = [1, 1, 3, 2, 1, 2]
  partition[:, i += 1] = [1, 1, 2, 2, 3, 1]
  partition[:, i += 1] = [1, 1, 2, 3, 2, 1]
  partition[:, i += 1] = [1, 1, 3, 2, 2, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 2, 3]
  partition[:, i += 1] = [1, 2, 1, 1, 3, 2]
  partition[:, i += 1] = [1, 3, 1, 1, 2, 2]
  partition[:, i += 1] = [1, 2, 1, 2, 1, 3]
  partition[:, i += 1] = [1, 2, 1, 3, 1, 2]
  partition[:, i += 1] = [1, 3, 1, 2, 1, 2]
  partition[:, i += 1] = [1, 2, 1, 2, 3, 1]
  partition[:, i += 1] = [1, 2, 1, 3, 2, 1]
  partition[:, i += 1] = [1, 3, 1, 2, 2, 1]
  partition[:, i += 1] = [1, 2, 2, 1, 1, 3]
  partition[:, i += 1] = [1, 2, 3, 1, 1, 2]
  partition[:, i += 1] = [1, 3, 2, 1, 1, 2]
  partition[:, i += 1] = [1, 2, 2, 1, 3, 1]
  partition[:, i += 1] = [1, 2, 3, 1, 2, 1]
  partition[:, i += 1] = [1, 3, 2, 1, 2, 1]
  partition[:, i += 1] = [1, 2, 2, 3, 1, 1]
  partition[:, i += 1] = [1, 2, 3, 2, 1, 1]
  partition[:, i += 1] = [1, 3, 2, 2, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1, 2, 3]
  partition[:, i += 1] = [2, 1, 1, 1, 3, 2]
  partition[:, i += 1] = [3, 1, 1, 1, 2, 2]
  partition[:, i += 1] = [2, 1, 1, 2, 1, 3]
  partition[:, i += 1] = [2, 1, 1, 3, 1, 2]
  partition[:, i += 1] = [3, 1, 1, 2, 1, 2]
  partition[:, i += 1] = [2, 1, 1, 2, 3, 1]
  partition[:, i += 1] = [2, 1, 1, 3, 2, 1]
  partition[:, i += 1] = [3, 1, 1, 2, 2, 1]
  partition[:, i += 1] = [2, 1, 2, 1, 1, 3]
  partition[:, i += 1] = [2, 1, 3, 1, 1, 2]
  partition[:, i += 1] = [3, 1, 2, 1, 1, 2]
  partition[:, i += 1] = [2, 1, 2, 1, 3, 1]
  partition[:, i += 1] = [2, 1, 3, 1, 2, 1]
  partition[:, i += 1] = [3, 1, 2, 1, 2, 1]
  partition[:, i += 1] = [2, 1, 2, 3, 1, 1]
  partition[:, i += 1] = [2, 1, 3, 2, 1, 1]
  partition[:, i += 1] = [3, 1, 2, 2, 1, 1]
  partition[:, i += 1] = [2, 2, 1, 1, 1, 3]
  partition[:, i += 1] = [2, 3, 1, 1, 1, 2]
  partition[:, i += 1] = [3, 2, 1, 1, 1, 2]
  partition[:, i += 1] = [2, 2, 1, 1, 3, 1]
  partition[:, i += 1] = [2, 3, 1, 1, 2, 1]
  partition[:, i += 1] = [3, 2, 1, 1, 2, 1]
  partition[:, i += 1] = [2, 2, 1, 3, 1, 1]
  partition[:, i += 1] = [2, 3, 1, 2, 1, 1]
  partition[:, i += 1] = [3, 2, 1, 2, 1, 1]
  partition[:, i += 1] = [2, 2, 3, 1, 1, 1]
  partition[:, i += 1] = [2, 3, 2, 1, 1, 1]
  partition[:, i += 1] = [3, 2, 2, 1, 1, 1]

  # **|**|**
  j += 1
  index[j] = i + 1
  k[j] = 3
  blocksize[:,  j] = [2, 2, 2, 0, 0, 0]
  blockcount[:, j] = [0, 3, 0, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 2, 2, 3, 3]
  partition[:, i += 1] = [1, 1, 2, 3, 2, 3]
  partition[:, i += 1] = [1, 1, 2, 3, 3, 2]
  partition[:, i += 1] = [1, 2, 1, 2, 3, 3]
  partition[:, i += 1] = [1, 2, 1, 3, 2, 3]
  partition[:, i += 1] = [1, 2, 1, 3, 3, 2]
  partition[:, i += 1] = [1, 2, 2, 1, 3, 3]
  partition[:, i += 1] = [1, 2, 3, 1, 2, 3]
  partition[:, i += 1] = [1, 2, 3, 1, 3, 2]
  partition[:, i += 1] = [1, 2, 2, 3, 1, 3]
  partition[:, i += 1] = [1, 2, 3, 2, 1, 3]
  partition[:, i += 1] = [1, 2, 3, 3, 1, 2]
  partition[:, i += 1] = [1, 2, 2, 3, 3, 1]
  partition[:, i += 1] = [1, 2, 3, 2, 3, 1]
  partition[:, i += 1] = [1, 2, 3, 3, 2, 1]

  # k = 4

  # ***|*|*|*
  j += 1
  index[j] = i + 1
  k[j] = 4
  blocksize[:,  j] = [3, 1, 1, 1, 0, 0]
  blockcount[:, j] = [3, 0, 1, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 1, 2, 3, 4]
  partition[:, i += 1] = [1, 1, 2, 1, 3, 4]
  partition[:, i += 1] = [1, 1, 2, 3, 1, 4]
  partition[:, i += 1] = [1, 1, 2, 3, 4, 1]
  partition[:, i += 1] = [1, 2, 1, 1, 3, 4]
  partition[:, i += 1] = [1, 2, 1, 3, 1, 4]
  partition[:, i += 1] = [1, 2, 1, 3, 4, 1]
  partition[:, i += 1] = [1, 2, 3, 1, 1, 4]
  partition[:, i += 1] = [1, 2, 3, 1, 4, 1]
  partition[:, i += 1] = [1, 2, 3, 4, 1, 1]
  partition[:, i += 1] = [2, 1, 1, 1, 3, 4]
  partition[:, i += 1] = [2, 1, 1, 3, 1, 4]
  partition[:, i += 1] = [2, 1, 1, 3, 4, 1]
  partition[:, i += 1] = [2, 1, 3, 1, 1, 4]
  partition[:, i += 1] = [2, 1, 3, 1, 4, 1]
  partition[:, i += 1] = [2, 1, 3, 4, 1, 1]
  partition[:, i += 1] = [2, 3, 1, 1, 1, 4]
  partition[:, i += 1] = [2, 3, 1, 1, 4, 1]
  partition[:, i += 1] = [2, 3, 1, 4, 1, 1]
  partition[:, i += 1] = [2, 3, 4, 1, 1, 1]

  # **|**|*|*
  j += 1
  index[j] = i + 1
  k[j] = 4
  blocksize[:,  j] = [2, 2, 1, 1, 0, 0]
  blockcount[:, j] = [2, 2, 0, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 2, 2, 3, 4]
  partition[:, i += 1] = [1, 1, 2, 3, 2, 4]
  partition[:, i += 1] = [1, 1, 2, 3, 4, 2]
  partition[:, i += 1] = [1, 1, 3, 2, 2, 4]
  partition[:, i += 1] = [1, 1, 3, 2, 4, 2]
  partition[:, i += 1] = [1, 1, 3, 4, 2, 2]
  partition[:, i += 1] = [1, 2, 1, 2, 3, 4]
  partition[:, i += 1] = [1, 2, 1, 3, 2, 4]
  partition[:, i += 1] = [1, 2, 1, 3, 4, 2]
  partition[:, i += 1] = [1, 3, 1, 2, 2, 4]
  partition[:, i += 1] = [1, 3, 1, 2, 4, 2]
  partition[:, i += 1] = [1, 3, 1, 4, 2, 2]
  partition[:, i += 1] = [1, 2, 2, 1, 3, 4]
  partition[:, i += 1] = [1, 2, 3, 1, 2, 4]
  partition[:, i += 1] = [1, 2, 3, 1, 4, 2]
  partition[:, i += 1] = [1, 3, 2, 1, 2, 4]
  partition[:, i += 1] = [1, 3, 2, 1, 4, 2]
  partition[:, i += 1] = [1, 3, 4, 1, 2, 2]
  partition[:, i += 1] = [1, 2, 2, 3, 1, 4]
  partition[:, i += 1] = [1, 2, 3, 2, 1, 4]
  partition[:, i += 1] = [1, 2, 3, 4, 1, 2]
  partition[:, i += 1] = [1, 3, 2, 2, 1, 4]
  partition[:, i += 1] = [1, 3, 2, 4, 1, 2]
  partition[:, i += 1] = [1, 3, 4, 2, 1, 2]
  partition[:, i += 1] = [1, 2, 2, 3, 4, 1]
  partition[:, i += 1] = [1, 2, 3, 2, 4, 1]
  partition[:, i += 1] = [1, 2, 3, 4, 2, 1]
  partition[:, i += 1] = [1, 3, 2, 2, 4, 1]
  partition[:, i += 1] = [1, 3, 2, 4, 2, 1]
  partition[:, i += 1] = [1, 3, 4, 2, 2, 1]
  partition[:, i += 1] = [3, 1, 1, 2, 2, 4]
  partition[:, i += 1] = [3, 1, 1, 2, 4, 2]
  partition[:, i += 1] = [3, 1, 1, 4, 2, 2]
  partition[:, i += 1] = [3, 1, 2, 1, 2, 4]
  partition[:, i += 1] = [3, 1, 2, 1, 4, 2]
  partition[:, i += 1] = [3, 1, 4, 1, 2, 2]
  partition[:, i += 1] = [3, 1, 2, 2, 1, 4]
  partition[:, i += 1] = [3, 1, 2, 4, 1, 2]
  partition[:, i += 1] = [3, 1, 4, 2, 1, 2]
  partition[:, i += 1] = [3, 1, 2, 2, 4, 1]
  partition[:, i += 1] = [3, 1, 2, 4, 2, 1]
  partition[:, i += 1] = [3, 1, 4, 2, 2, 1]
  partition[:, i += 1] = [3, 4, 1, 1, 2, 2]
  partition[:, i += 1] = [3, 4, 1, 2, 1, 2]
  partition[:, i += 1] = [3, 4, 1, 2, 2, 1]

  # k = 5

  # **|*|*|*|*
  j += 1
  index[j] = i + 1
  k[j] = 5
  blocksize[:,  j] = [2, 1, 1, 1, 1, 0]
  blockcount[:, j] = [4, 1, 0, 0, 0, 0]

  partition[:, i += 1] = [1, 1, 2, 3, 4, 5]
  partition[:, i += 1] = [1, 2, 1, 3, 4, 5]
  partition[:, i += 1] = [1, 2, 3, 1, 4, 5]
  partition[:, i += 1] = [1, 2, 3, 4, 1, 5]
  partition[:, i += 1] = [1, 2, 3, 4, 5, 1]
  partition[:, i += 1] = [2, 1, 1, 3, 4, 5]
  partition[:, i += 1] = [2, 1, 3, 1, 4, 5]
  partition[:, i += 1] = [2, 1, 3, 4, 1, 5]
  partition[:, i += 1] = [2, 1, 3, 4, 5, 1]
  partition[:, i += 1] = [2, 3, 1, 1, 4, 5]
  partition[:, i += 1] = [2, 3, 1, 4, 1, 5]
  partition[:, i += 1] = [2, 3, 1, 4, 5, 1]
  partition[:, i += 1] = [2, 3, 4, 1, 1, 5]
  partition[:, i += 1] = [2, 3, 4, 1, 5, 1]
  partition[:, i += 1] = [2, 3, 4, 5, 1, 1]

  # k = 6

  # *|*|*|*|*|*
  j += 1
  index[j] = i + 1
  k[j] = 6
  blocksize[:,  j] = [1, 1, 1, 1, 1, 1]
  blockcount[:, j] = [6, 0, 0, 0, 0, 0]

  partition[:, i += 1] = [1, 2, 3, 4, 5, 6]

  TestPartition(n, B, C, k, index, blocksize, blockcount, partition)
end

function TestPartition(n::Int)
  if n == 1
    TestP1()
  elseif n == 2
    TestP2()
  elseif n == 3
    TestP3()
  elseif n == 4
    TestP4()
  elseif n == 5
    TestP5()
  elseif n == 6
    TestP6()
  else
    throw(Kpax3DomainError("Argument n must be positive and lesser than 7."))
  end
end
