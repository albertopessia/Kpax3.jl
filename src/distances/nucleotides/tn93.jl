# This file is part of Kpax3. License is MIT.

#=
References:

Tamura K. and Nei M. (1993). Estimation of the number of nucleotide
substitutions in the control region of mitochondrial DNA in humans and
chimpanzees. Mol Biol Evol 10 (3): 512-526.
http://mbe.oxfordjournals.org/content/10/3/512

Tamura K. and Kumar S. (2002). Evolutionary Distance Estimation Under
Heterogeneous Substitution Pattern Among Lineages. Mol Biol Evol 19 (10):
1727-1736. http://mbe.oxfordjournals.org/content/19/10/1727
=#

#=
Description:

Compute Tamura Nei (1993) pairwise distances of dna sequences.

Arguments:

  data::Matrix{UInt8}
    m-by-n data matrix, where m is the common sequence length and n is the
    sample size
  ref::Vector{UInt8}
    reference sequence, i.e. a vector of length m storing the values of
    homogeneous sites

Details:

Only the four basic nucleotides are considered in the computations. It is
expected that Uracil has a value of 4 (equal to Thymidine).

If a pairwise distance is equal to -1.0, it means that is wasn't possible to
compute it. This usually happens when the hypotheses of the underlying
evolutionary model are not satisfied.

Value:

  d::Vector{Float64}
    evolutionary distances. Vector length is equal to n * (n - 1) / 2. It
    contains the values of the lower triangular matrix of the full distance
    matrix, ordered by column.
    To access the distance between units i and j (i < j), use
    d[n * (i - 1) - div(i * (i - 1), 2) + j - i]
=#
function distnttn93(data::Matrix{UInt8},
                    ref::Vector{UInt8})
  n = size(data, 2)

  d = zeros(Float64, div(n * (n - 1), 2))

  gt = zeros(Float64, 4)
  gb = zeros(Float64, 4)

  for b in 1:m
    if ref[b] == 1
      gt[1] += 1
    elseif ref[b] == 2
      gt[2] += 1
    elseif ref[b] == 3
      gt[3] += 1
    elseif ref[b] == 4
      gt[4] += 1
    end
  end

  for a in 1:n, b in 1:m
    if data[b, a] == 1
      gb[1] += 1
    elseif data[b, a] == 2
      gb[2] += 1
    elseif data[b, a] == 3
      gb[3] += 1
    elseif data[b, a] == 4
      gb[4] += 1
    end
  end

  gb[1] += n * gt[1]
  gb[2] += n * gt[2]
  gb[3] += n * gt[3]
  gb[4] += n * gt[4]

  tot = gb[1] + gb[2] + gb[3] + gb[4]

  gb[1] /= tot
  gb[2] /= tot
  gb[3] /= tot
  gb[4] /= tot

  h = gt[1] + gt[2] + gt[3] + gt[4]

  gr = gb[1] + gb[3]
  gy = gb[2] + gb[4]

  k1 = (gb[1] * gb[3]) / gr
  k2 = (gb[2] * gb[4]) / gy
  k3 = gr * gy - k1 * gy - k2 * gr

  idx = 0
  for j in 1:(n - 1), i in (j + 1):n
    d[idx += 1] = nttn93(i, j, data, gr, gy, h, k1, k2, k3)
  end

  d
end

#=
Description:

Compute Tamura Nei (1993) pairwise distances of dna sequences, with the Tamura
and Kumar (2002) correction for heterogeneous patterns.

Arguments:

  data::Matrix{UInt8}
    m-by-n data matrix, where m is the common sequence length and n is the
    sample size
  ref::Vector{UInt8}
    reference sequence, i.e. a vector of length m storing the values of
    homogeneous sites

Details:

Only the four basic nucleotides are considered in the computations. It is
expected that Uracil has a value of 4 (equal to Thymidine).

If a pairwise distance is equal to -1.0, it means that is wasn't possible to
compute it. This usually happens when the hypotheses of the underlying
evolutionary model are not satisfied.

Value:

  d::Vector{Float64}
    evolutionary distances. Vector length is equal to n * (n - 1) / 2. It
    contains the values of the lower triangular matrix of the full distance
    matrix, ordered by column.
    To access the distance between units i and j (i < j), use
    d[n * (i - 1) - div(i * (i - 1), 2) + j - i]
=#
function distntmtn93(data::Matrix{UInt8},
                     ref::Vector{UInt8})
  n = size(data, 2)

  d = zeros(Float64, div(n * (n - 1), 2))

  gt = zeros(Float64, 4)
  gb = zeros(Float64, 4)

  for b in 1:m
    if ref[b] == 1
      gt[1] += 1
    elseif ref[b] == 2
      gt[2] += 1
    elseif ref[b] == 3
      gt[3] += 1
    elseif ref[b] == 4
      gt[4] += 1
    end
  end

  for a in 1:n, b in 1:m
    if data[b, a] == 1
      gb[1] += 1
    elseif data[b, a] == 2
      gb[2] += 1
    elseif data[b, a] == 3
      gb[3] += 1
    elseif data[b, a] == 4
      gb[4] += 1
    end
  end

  gb[1] += n * gt[1]
  gb[2] += n * gt[2]
  gb[3] += n * gt[3]
  gb[4] += n * gt[4]

  tot = gb[1] + gb[2] + gb[3] + gb[4]

  gb[1] /= tot
  gb[2] /= tot
  gb[3] /= tot
  gb[4] /= tot

  h = gt[1] + gt[2] + gt[3] + gt[4]

  gr = gb[1] + gb[3]
  gy = gb[2] + gb[4]

  k1 = (gb[1] * gb[3]) / gr
  k2 = (gb[2] * gb[4]) / gy
  k3 = gr * gy - k1 * gy - k2 * gr

  idx = 0
  for j in 1:(n - 1), i in (j + 1):n
    d[idx += 1] = ntmtn93(i, j, data, gt, gr, gy, h, k1, k2, k3)
  end

  d
end

#=
Description:

Basic function used to compute the Tamura Nei (1993) distance between two dna
sequences.

Arguments:

  i::Int
    index of first sequence
  j::Int
    index of second sequence
  data::Matrix{UInt8}
    m-by-n data matrix, where m is the common sequence length and n is the
    sample size
  gr::Float64
    gr = pA + pG. Observed proportion (in the whole dataset) of purines
  gy::Float64
    gy = pC + pT. Observed proportion (in the whole dataset) of pyrimidines
  h::Float64
    total number of homogeneous sites that are not missing
  k1::Float64
    k1 = (pA * pG) / gr. pA and pG are the proportions of nucleotide A and
    nucleotide G observed in the whole dataset, respectively
  k2::Float64
    k2 = (pC * pT) / gy. pC and pT are the proportions of nucleotide C and
    nucleotide T observed in the whole dataset, respectively
  k3::Float64
    k3 = gr * gy - k1 * gy - k2 * gr

Value:

  d::Float64
    evolutionary distance between the two dna sequences
=#
function nttn93(i::Int,
                j::Int,
                data::Matrix{UInt8},
                gr::Float64,
                gy::Float64,
                h::Float64,
                k1::Float64,
                k2::Float64,
                k3::Float64)
  d = -1.0

  # proportions of transitional differences between nucleotides A and G
  ag = 0.0

  # proportions of transitional differences between nucleotides C and T
  ct = 0.0

  # proportions of transversional differences
  ry = 0.0

  x1 = 0x00
  x2 = 0x00

  for b in 1:size(data, 1)
    x1 = data[b, i]
    x2 = data[b, j]

    if (0x00 < x1 < 0x05) && (0x00 < x2 < 0x05)
      if ((x1 == 0x01) && (x2 == 0x03)) || ((x1 == 0x03) && (x2 == 0x01))
        ag += 1
      elseif ((x1 == 0x02) && (x2 == 0x04)) || ((x1 == 0x04) && (x2 == 0x02))
        ct += 1
      elseif (((x1 == 0x01) || (x1 == 0x03)) &&
              ((x2 == 0x02) || (x2 == 0x04))) ||
             (((x1 == 0x02) || (x1 == 0x04)) &&
              ((x2 == 0x01) || (x2 == 0x03)))
        ry += 1
      end

      h += 1
    end
  end

  if (ag + ct + ry) > 0
    ag /= h
    ct /= h
    ry /= h

    w1 = 1 - (ag / k1 + ry / gr) / 2
    w2 = 1 - (ct / k2 + ry / gy) / 2
    w3 = 1 - ry / (2 * gr * gy)

    if (w1 > 0) && (w2 > 0) && (w3 > 0)
      d = - 2 * (k1 * log(w1) + k2 * log(w2) + k3 * log(w3))
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

Basic function used to compute the Tamura Nei (1993) distance between two dna
sequences, with the Tamura and Kumar (2002) correction for heterogeneous
patterns.

Arguments:

  i::Int
    index of first sequence
  j::Int
    index of second sequence
  data::Matrix{UInt8}
    m-by-n data matrix, where m is the common sequence length and n is the
    sample size
  gt::Vector{Float64}
    common absolute frequency for the 4 nucleotides, i.e. the count of each
    nucleotide at homogeneous sites
  gr::Float64
    gr = pA + pG. Observed proportion (in the whole dataset) of purines
  gy::Float64
    gy = pC + pT. Observed proportion (in the whole dataset) of pyrimidines
  h::Float64
    total number of homogeneous sites that are not missing
  k1::Float64
    k1 = (pA * pG) / gr. pA and pG are the proportions of nucleotide A and
    nucleotide G observed in the whole dataset, respectively
  k2::Float64
    k2 = (pC * pT) / gy. pC and pT are the proportions of nucleotide C and
    nucleotide T observed in the whole dataset, respectively
  k3::Float64
    k3 = gr * gy - k1 * gy - k2 * gr

Value:

  d::Float64
    evolutionary distance between the two dna sequences
=#
function ntmtn93(i::Int,
                 j::Int,
                 data::Matrix{UInt8},
                 gt::Vector{Float64},
                 gr::Float64,
                 gy::Float64,
                 h::Float64,
                 k1::Float64,
                 k2::Float64,
                 k3::Float64)
  d = -1.0

  # effective length, i.e. total number of sites at which both sequences have
  # non-missing values
  n = fill(h, 3)

  # proportion of observed nucleotides
  g1 = copy(gt)
  g2 = copy(gt)

  # proportions of transitional differences between nucleotides A and G
  ag = 0.0

  # proportions of transitional differences between nucleotides C and T
  ct = 0.0

  # proportions of transversional differences
  ry = 0.0

  r1 = 0.0
  y1 = 0.0

  r2 = 0.0
  y2 = 0.0

  x1 = 0x00
  x2 = 0x00

  for b in 1:size(data, 1)
    x1 = data[b, i]
    x2 = data[b, j]

    if 0x00 < x1 < 0x05
      g1[x1] += 1
      n[1] += 1
    end

    if 0x00 < x2 < 0x05
      g2[x2] += 1
      n[2] += 1
    end

    if (0x00 < x1 < 0x05) && (0x00 < x2 < 0x05)
      if ((x1 == 0x01) && (x2 == 0x03)) || ((x1 == 0x03) && (x2 == 0x01))
        ag += 1
      elseif ((x1 == 0x02) && (x2 == 0x04)) || ((x1 == 0x04) && (x2 == 0x02))
        ct += 1
      elseif (((x1 == 0x01) || (x1 == 0x03)) &&
              ((x2 == 0x02) || (x2 == 0x04))) ||
             (((x1 == 0x02) || (x1 == 0x04)) &&
              ((x2 == 0x01) || (x2 == 0x03)))
        ry += 1
      end

      n[3] += 1
    end
  end

  if (ag + ct + ry) > 0
    g1 /= n[1]
    g2 /= n[2]

    ag /= n[3]
    ct /= n[3]
    ry /= n[3]

    r1 = g1[1] + g1[3]
    y1 = g1[2] + g1[4]

    r2 = g2[1] + g2[3]
    y2 = g2[2] + g2[4]

    w1 = 1 - (gr * ag) / (g1[1] * g2[3] + g1[3] * g2[1]) - ry / (2 * gr)
    w2 = 1 - (gy * ct) / (g1[2] * g2[4] + g1[4] * g2[2]) - ry / (2 * gy)
    w3 = 1 - ry / (r1 * y2 + y1 * r2)

    if (w1 > 0) && (w2 > 0) && (w3 > 0)
      d = - 2 * (k1 * log(w1) + k2 * log(w2) + k3 * log(w3))
    end
  else
    # sequences are identical (or couldn't be compared because of missing
    # values. But in this case, by default, we consider them identical)
    d = 0.0
  end

  d
end
