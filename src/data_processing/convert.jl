# This file is part of Kpax3. License is MIT.

"""
# Convert categorical (integer) data to binary

## Description

Convert an integer matrix to a binary (indicator) matrix.

## Usage

categorical2binary(data, maxval)

## Arguments

* `data` Integer matrix to be converted
* `maxval` Theoretical maximum value observable in `data`

## Value

A tuple containing the following variables:

* `bindata` Original data matrix encoded as a binary (indicator) matrix
* `val` vector with unique values per MSA site
* `key` vector with indices of each value

## Example

If `data` consists of just the following three units

    0 2 1
    1 3 2
    2 4 0
    2 2 3

then `bindata` will be equal to

    0 0 1
    0 1 0
    1 0 0
    0 0 1
    0 1 0
    1 0 0
    0 1 0
    1 1 0
    0 0 1

while

    val = [1, 2, 1, 2, 3, 2, 4, 2, 3] (i.e. 12 123 24 23)
    key = [1, 1, 2, 2, 2, 3, 3, 4, 4] (i.e. 11 222 33 44)

`0` values (missing data) are discarded.
"""
function categorical2binary{T <: Integer}(data::Matrix{T},
                                          maxval::T)
  # TOOPTIMIZE: traversing the matrices by row instead of by column
  n = size(data, 2)
  m = 0

  # number of unique values per column
  v = zeros(Int, size(data, 1))

  # maximum value actually observed in data
  M = c = zero(T)

  tmp = falses(Int(maxval))
  for row in 1:size(data, 1)
    fill!(tmp, false)

    for col in 1:n
      c = data[row, col]

      if (c < 0) || (c > maxval)
        throw(KDomainError(string("Value outside the allowed range  at (row, ",
                                  "col) = (", row, ", ", col, ").")))
      end

      if (c > 0) && !tmp[c]
        tmp[c] = true
        v[row] += 1

        if c > M
          M = c
        end
      end
    end

    m += v[row]
  end

  bindata = zeros(UInt8, m, n)
  val = zeros(T, m)
  key = zeros(Int, m)

  i = 0
  tmp1 = falses(Int(M))
  tmp2 = zeros(UInt8, M, n)
  for row in 1:size(data, 1)
    fill!(tmp1, false)
    fill!(tmp2, false)

    for col in 1:n
      c = data[row, col]

      if c > 0
        tmp1[c] = true
        tmp2[c, col] = 1
      end
    end

    bindata[(i + 1):(i + v[row]), :] = tmp2[tmp1, :]
    val[(i + 1):(i + v[row])] = find(tmp1)
    key[(i + 1):(i + v[row])] = row

    i += v[row]
  end

  (bindata, val, key)
end
