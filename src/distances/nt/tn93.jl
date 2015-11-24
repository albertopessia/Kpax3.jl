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
function distnttn93(rawdata::Array{Uint8, 2},
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

  gr = gb[1] + gb[3]
  gy = gb[2] + gb[4]

  k1 = (gb[1] * gb[3]) / gr
  k2 = (gb[2] * gb[4]) / gy
  k3 = gr * gy - k1 * gy - k2 * gr

  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    d[idx] = nttn93(rawdata[:, i], rawdata[:, j], h, gr, gy, k1, k2, k3)
    idx += 1
  end

  d
end

#=
Description:

Compute Tamura Nei (1993) pairwise distances of dna sequences, with the Tamura
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
function distntmtn93(rawdata::Array{Uint8, 2},
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

  gr = gb[1] + gb[3]
  gy = gb[2] + gb[4]

  k1 = (gb[1] * gb[3]) / gr
  k2 = (gb[2] * gb[4]) / gy
  k3 = gr * gy - k1 * gy - k2 * gr

  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    d[idx] = ntmtn93(rawdata[:, i], rawdata[:, j], gt, h, gr, gy, k1, k2, k3)
    idx += 1
  end

  d
end

#=
Description:

Basic function used to compute the Tamura Nei (1993) distance between two dna
sequences.

Arguments:

  s1::Array{Uint8, 1}
    dna sequence
  s2::Array{Uint8, 1}
    dna sequence
  n::Float64
    total number of non-missing homogeneous sites
  gr::Float64
    gr = pA + pG. Observed proportion (in the whole dataset) of purines
  gy::Float64
    gy = pC + pT. Observed proportion (in the whole dataset) of pyrimidines
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
function nttn93(x1::Array{Uint8, 1},
                x2::Array{Uint8, 1},
                n::Float64,
                gr::Float64,
                gy::Float64,
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

  for i in 1:length(x1)
    if (0 < x1[i] < 5) && (0 < x2[i] < 5)
      if ((x1[i] == 1) && (x2[i] == 3)) || ((x1[i] == 3) && (x2[i] == 1))
        ag += 1
      elseif ((x1[i] == 2) && (x2[i] == 4)) || ((x1[i] == 4) && (x2[i] == 2))
        ct += 1
      elseif (((x1[i] == 1) || (x1[i] == 3)) &&
              ((x2[i] == 2) || (x2[i] == 4))) ||
             (((x1[i] == 2) || (x1[i] == 4)) &&
              ((x2[i] == 1) || (x2[i] == 3)))
        ry += 1
      end

      n += 1
    end
  end

  if (ag + ct + ry) > 0
    ag /= n
    ct /= n
    ry /= n

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
  gr::Float64
    gr = pA + pG. Observed proportion (in the whole dataset) of purines
  gy::Float64
    gy = pC + pT. Observed proportion (in the whole dataset) of pyrimidines
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
function ntmtn93(x1::Array{Uint8, 1},
                 x2::Array{Uint8, 1},
                 gt::Array{Float64, 1},
                 h::Float64,
                 gr::Float64,
                 gy::Float64,
                 k1::Float64,
                 k2::Float64,
                 k3::Float64)
  d = -1.0

  # effective length, i.e. total number of sites at which both sequences have
  # non-missing values
  n = zeros(Float64, 3)
  n[:] = h

  # proportion of observed nucleotides
  g1 = zeros(Float64, 4)
  g1[:] = gt

  g2 = zeros(Float64, 4)
  g2[:] = gt

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
      if ((x1[i] == 1) && (x2[i] == 3)) || ((x1[i] == 3) && (x2[i] == 1))
        ag += 1
      elseif ((x1[i] == 2) && (x2[i] == 4)) || ((x1[i] == 4) && (x2[i] == 2))
        ct += 1
      elseif (((x1[i] == 1) || (x1[i] == 3)) &&
              ((x2[i] == 2) || (x2[i] == 4))) ||
             (((x1[i] == 2) || (x1[i] == 4)) &&
              ((x2[i] == 1) || (x2[i] == 3)))
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
