# This file is part of Kpax3. License is MIT.

"""
# Convert categorical (integer) data to binary

## Description

Convert an integer matrix to a binary (indicator) matrix.

## Usage

categorical2binary(data)

## Arguments

* `data` Integer matrix to be converted

## Value

A tuple containing the following variables:

* `bindata` Original data matrix encoded as a binary (indicator) matrix
* `val` vector with unique values per MSA site
* `key` vector with indices of each value

## Example

If `data` consists of just the following three rows

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

Missing data (`0` values) is discarded.
"""
function categorical2binary{T <: Integer}(data::Array{T, 2})
  # TOOPTIMIZE: traversing the matrices by row instead of by column
  n = size(data, 2)
  m = 0

  # number of unique values per column
  v = zeros(Int, size(data, 1))

  # maximum value observed in data
  maxval = c = zero(T)

  tmp1 = falses(n)
  for row in 1:size(data, 1)
    tmp1[:] = false

    for col in 1:size(data, 2)
      c = data[row, col]

      if (c > 0) && !tmp1[c]
        tmp1[c] = true
        v[row] += 1

        if c > maxval
          maxval = c
        end
      end
    end

    m += v[row]
  end

  bindata = zeros(UInt8, m, n)
  val = zeros(Int, m)
  key = zeros(Int, m)

  i = 0
  tmp1 = falses(Int(maxval))
  tmp2 = zeros(UInt8, maxval, n)
  for row in 1:size(data, 1)
    tmp1[:] = false
    tmp2[:] = false

    for col in 1:size(data, 2)
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
