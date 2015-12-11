# This file is part of Kpax3. License is MIT.

#=
References:

Tajima F. and Nei M. (1984). Estimation of evolutionary distance between
nucleotide sequences. Mol Biol Evol 1 (3):269-285.
http://mbe.oxfordjournals.org/content/1/3/269

Tamura K. and Kumar S. (2002). Evolutionary Distance Estimation Under
Heterogeneous Substitution Pattern Among Lineages. Mol Biol Evol 19 (10):
1727-1736. http://mbe.oxfordjournals.org/content/19/10/1727
=#

#=
Description:

Compute Tajima Nei (1984) pairwise distances of protein sequences, with the
Tamura and Kumar (2002) correction for heterogeneous patterns.

Arguments:

  rawdata::Array{UInt8, 2}
    m-by-n data matrix, where m is the common sequence length and n is the
    sample size
  ref::Array{UInt8, 1}
    reference sequence, i.e. a vector of length m storing the values of
    homogeneous sites

Details:

Pyrrolysine (O) and Selenocysteine (U) are treated as valid amino acids if they
have not been given a value of zero (missing).

If a pairwise distance is equal to -1.0, it means that is wasn't possible to
compute it. This usually happens when the hypotheses of the underlying
evolutionary model are not satisfied.

Value:

  d::Array{Float64, 1}
    evolutionary distances. Vector length is equal to n * (n - 1) / 2. It
    contains the values of the lower triangular matrix of the full distance
    matrix, ordered by column.
    To access the distance between units i and j (i < j), use
    d[n * (i - 1) - div(i * (i - 1), 2) + j - i]
=#
function distaamtn84(data::Matrix{UInt8},
                     ref::Vector{UInt8})
  m, n = size(data)

  d = zeros(Float64, div(n * (n - 1), 2))

  tmpref = zeros(Float64, 30)
  tmpraw = zeros(Float64, 29)

  for i in 1:length(ref)
    tmpref[ref[i] + 1] += 1
  end

  for col in 1:n, row in 1:m
    tmpraw[data[row, col] + 1] += 1
  end

  gt = tmpref[2:23]
  gb = tmpraw[2:23] + n * gt

  st = 0.0
  sb = 0.0
  for aa in 1:22
    st += gt[aa]
    sb += gb[aa]
  end

  gb /= sb

  sc = 0.0
  for aa in 1:22
    sc += gb[aa]^2
  end

  b = 1 - sc

  idx = 0
  for j in 1:(n - 1), i in (j + 1):n
    d[idx += 1] = aamtn84(i, j, data, gt, st, b)
  end

  d
end

#=
Description:

Compute the Tajima Nei (1984) distance between two protein sequences, with the
Tamura and Kumar (2002) correction for heterogeneous patterns.

Arguments:

  s1::Array{UInt8, 1}
    protein sequence
  s2::Array{UInt8, 1}
    protein sequence
  gt::Array{Float64, 1}
    common absolute frequency for the 22 amino acids, i.e. the count of each
    amino acid at homogeneous sites
  h::Float64
    h = sum(gt), i.e. the total number of homogeneous sites that are not missing
  b::Float64
    b = 1 - sum(pi^2), where pi is the proportion of amino acid i observed in
    the whole dataset

Value:

  d::Float64
    evolutionary distance between the two protein sequences
=#
function aamtn84(i::Int,
                 j::Int,
                 data::Matrix{UInt8},
                 gt::Vector{Float64},
                 h::Float64,
                 b::Float64)
  d = -1.0

  # effective length, i.e. total number of sites at which both sequences have
  # non-missing values
  n = fill(h, 3)

  # proportion of different elements
  p = 0.0

  # proportion of observed amino acids
  g1 = copy(gt)
  g2 = copy(gt)

  f = 0.0
  w = 0.0

  x1 = 0x00
  x2 = 0x00

  for b in 1:size(data, 1)
    x1 = data[b, i]
    x2 = data[b, j]

    if UInt8(0) < x1 < UInt8(23)
      g1[x1] += 1
      n[1] += 1
    end

    if UInt8(0) < x2 < UInt8(23)
      g2[x2] += 1
      n[2] += 1
    end

    if (UInt8(0) < x1 < UInt8(23)) && (UInt8(0) < x2 < UInt8(23))
      if x1 != x2
        p += 1
      end

      n[3] += 1
    end
  end

  if p > 0
    g1 /= n[1]
    g2 /= n[2]

    p /= n[3]

    for i in 1:22
      f += g1[i] * g2[i]
    end

    x = 1 - (p / (1 - f))

    if x > 0
      d = - b * log(x)
    end
  else
    # sequences are identical (or couldn't be compared because of missing
    # values. But in this case, by default, we consider them identical)
    d = 0.0
  end

  d
end
