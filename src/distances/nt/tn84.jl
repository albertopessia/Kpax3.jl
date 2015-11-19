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

Compute Tajima Nei (1984) pairwise distances of dna sequences.

Arguments:

  rawdata::Array{Uint8, 2}
    m-by-n data matrix, where m is the common sequence length and n is the
    sample size
  ref::Array{Uint8, 1}
    reference sequence, i.e. a vector of length m storing the values of
    homogeneous sites

Details:

Only the four basic nucleotides are considered in the computations. It is
expected that Uracil has a value of 4 (equal to Thymidine).

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
function distnttn84(rawdata::Array{Uint8, 2},
                    ref::Array{Uint8, 1})
  n = size(rawdata, 2)
  d = zeros(Float64, div(n * (n - 1), 2))

  gt = zeros(Float64, 4)
  gb = zeros(Float64, 4)
  gp = zeros(Float64, 6)

  gt[1] = sum(ref .== 1)
  gt[2] = sum(ref .== 2)
  gt[3] = sum(ref .== 3)
  gt[4] = sum(ref .== 4)

  gb[1] = sum(rawdata .== 1) + n * gt[1]
  gb[2] = sum(rawdata .== 2) + n * gt[2]
  gb[3] = sum(rawdata .== 3) + n * gt[3]
  gb[4] = sum(rawdata .== 4) + n * gt[4]

  gb /= (gb[1] + gb[2] + gb[3] + gb[4])

  h = gt[1] + gt[2] + gt[3] + gt[4]

  gp[1] = 2 * gb[1] * gb[2]
  gp[2] = 2 * gb[1] * gb[3]
  gp[3] = 2 * gb[1] * gb[4]
  gp[4] = 2 * gb[2] * gb[3]
  gp[5] = 2 * gb[2] * gb[4]
  gp[6] = 2 * gb[3] * gb[4]

  b1 = 1 - (gb[1]^2 + gb[2]^2 + gb[3]^2 + gb[4]^2)

  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    d[idx] = nttn84(rawdata[:, i], rawdata[:, j], h, gp, b1)
    idx += 1
  end

  d
end

#=
Description:

Compute Tajima Nei (1984) pairwise distances of dna sequences, with the Tamura
and Kumar (2002) correction for heterogeneous patterns.

Arguments:

  rawdata::Array{Uint8, 2}
    m-by-n data matrix, where m is the common sequence length and n is the
    sample size
  ref::Array{Uint8, 1}
    reference sequence, i.e. a vector of length m storing the values of
    homogeneous sites

Details:

Only the four basic nucleotides are considered in the computations. It is
expected that Uracil has a value of 4 (equal to Thymidine).

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
function distntmtn84(rawdata::Array{Uint8, 2},
                     ref::Array{Uint8, 1})
  n = size(rawdata, 2)
  d = zeros(Float64, div(n * (n - 1), 2))

  gt = zeros(Float64, 4)
  gb = zeros(Float64, 4)

  gt[1] = sum(ref .== 1)
  gt[2] = sum(ref .== 2)
  gt[3] = sum(ref .== 3)
  gt[4] = sum(ref .== 4)

  gb[1] = sum(rawdata .== 1) + n * gt[1]
  gb[2] = sum(rawdata .== 2) + n * gt[2]
  gb[3] = sum(rawdata .== 3) + n * gt[3]
  gb[4] = sum(rawdata .== 4) + n * gt[4]

  gb /= (gb[1] + gb[2] + gb[3] + gb[4])

  h = gt[1] + gt[2] + gt[3] + gt[4]

  b = 1 - (gb[1]^2 + gb[2]^2 + gb[3]^2 + gb[4]^2)

  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    d[idx] = ntmtn84(rawdata[:, i], rawdata[:, j], gt, h, b)
    idx += 1
  end

  d
end

#=
Description:

Compute the Tajima Nei (1984) distance between two dna sequences.

Arguments:

  s1::Array{Uint8, 1}
    dna sequence
  s2::Array{Uint8, 1}
    dna sequence
  n::Float64
    total number of non-missing homogeneous sites
  gp::Array{Float64, 1}
    gp_ij = 2 * pi * pj, where pi and pj are the proportions of nucleotides i
    and j, respectively, observed in the whole dataset
  b::Float64
    b = 1 - (pA^2 + pC^2 + pG^2 + pT^2), where pi is the proportion of
    nucleotide i observed in the whole dataset

Value:

  d::Float64
    evolutionary distance between the two dna sequences
=#
function nttn84(x1::Array{Uint8, 1},
                x2::Array{Uint8, 1},
                n::Float64,
                gp::Array{Float64, 1},
                b1::Float64)
  d = -1.0

  # proportion of different elements
  p = 0.0

  # proportion of nucleotide pairs
  pn = zeros(Float64, 6)

  c = 0.0
  w = 0.0

  for i in 1:length(x1)
    if (0 < x1[i] < 5) && (0 < x2[i] < 5)
      if x1[i] != x2[i]
        p += 1

        if     ((x1[i] == 1) && (x2[i] == 2)) || ((x1[i] == 2) && (x2[i] == 1))
          pn[1] += 1
        elseif ((x1[i] == 1) && (x2[i] == 3)) || ((x1[i] == 3) && (x2[i] == 1))
          pn[2] += 1
        elseif ((x1[i] == 1) && (x2[i] == 4)) || ((x1[i] == 4) && (x2[i] == 1))
          pn[3] += 1
        elseif ((x1[i] == 2) && (x2[i] == 3)) || ((x1[i] == 3) && (x2[i] == 2))
          pn[4] += 1
        elseif ((x1[i] == 2) && (x2[i] == 4)) || ((x1[i] == 4) && (x2[i] == 2))
          pn[5] += 1
        elseif ((x1[i] == 3) && (x2[i] == 4)) || ((x1[i] == 4) && (x2[i] == 3))
          pn[6] += 1
        end
      end

      n += 1
    end
  end

  if p > 0
    p /= n
    pn /= n

    c = pn[1]^2 / gp[1] + pn[2]^2 / gp[2] + pn[3]^2 / gp[3] +
        pn[4]^2 / gp[4] + pn[5]^2 / gp[5] + pn[6]^2 / gp[6]

    b = (b1 + p^2 / c) / 2

    w = 1 - p / b

    if w > 0
      d = - b * log(w)
    end
  else
    # sequences are identical (or couldn't be compared because of missing
    # values. But in this case, by default, we consider them identical)
    d = 0.0
  end

  d
end

#=
Description:

Compute the Tajima Nei (1984) distance between two dna sequences, with the
Tamura and Kumar (2002) correction for heterogeneous patterns.

Arguments:

  s1::Array{Uint8, 1}
    dna sequence
  s2::Array{Uint8, 1}
    dna sequence
  gt::Array{Float64, 1}
    common absolute frequency for the 4 nucleotides, i.e. the count of each
    nucleotide at homogeneous sites
  h::Float64
    h = gt[1] + gt[2] + gt[3] + gt[4], i.e. the total number of homogeneous
    sites that are not missing
  b::Float64
    b = 1 - (pA^2 + pC^2 + pG^2 + pT^2), where pi is the proportion of
    nucleotide i observed in the whole dataset

Value:

  d::Float64
    evolutionary distance between the two dna sequences
=#
function ntmtn84(x1::Array{Uint8, 1},
                 x2::Array{Uint8, 1},
                 gt::Array{Float64, 1},
                 h::Float64,
                 b::Float64)
  d = -1.0

  # effective length, i.e. total number of sites at which both sequences have
  # non-missing values
  n = zeros(Float64, 3)
  n[:] = h

  # proportion of different elements
  p = 0.0

  # proportion of observed nucleotides
  g1 = zeros(Float64, 4)
  g1[1:4] = gt

  g2 = zeros(Float64, 4)
  g2[1:4] = gt

  f = 0.0
  w = 0.0

  for i in 1:length(x1)
    if 0 < x1[i] < 5
      g1[x1[i]] += 1
      n[1] += 1
    end

    if 0 < x2[i] < 5
      g2[x2[i]] += 1
      n[2] += 1
    end

    if (0 < x1[i] < 5) && (0 < x2[i] < 5)
      if x1[i] != x2[i]
        p += 1
      end

      n[3] += 1
    end
  end

  if p > 0
    g1 /= n[1]
    g2 /= n[2]

    p /= n[3]

    f = 1 - (g1[1] * g2[1] + g1[2] * g2[2] + g1[3] * g2[3] + g1[4] * g2[4])
    w = 1 - p / f

    if w > 0
      d = - b * log(w)
    end
  else
    # sequences are identical (or couldn't be compared because of missing
    # values. But in this case, by default, we consider them identical)
    d = 0.0
  end

  d
end
