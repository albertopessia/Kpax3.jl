# This file is part of Kpax3. License is MIT.

fastafile = "data/mcmc.fasta"
partition = "data/mcmc.csv"
outfile = "../build/output.dat"
T = 1000000
burnin = 100000
tstep = 1
op = [1.0; 0.0; 0.0]
α = 0.0
θ = 1.0
γ = [0.6; 0.35; 0.05]
r = log(0.001) / log(0.95)
λs1 = 1.0
λs2 = 1.0
parawm = 5.0
maxclust = 1
maxunit = 1
verbose = true
verbosestep = 10000

# true parameters are
# R = [1; 1; 1; 1; 2; 2]
# S = [1; 1; 3; 3]
# C = [1 1 4 3;
#      1 1 3 4]

#=
# compute log(normalization constant) as accurate as possible
# first of all, find the posterior mode
include("data/partitions.jl")

function computelognormconst(cs,
                             k::Int,
                             lumpp::Float64,
                             data::Array{UInt8, 2},
                             po::TestPartition,
                             γ::Array{Float64, 1},
                             r::Float64,
                             priorR::PriorRowPartition)
  m, n = size(data)
  priorC = AminoAcidPriorCol(data, k, γ, r)

  st = po.index[po.k .== k][1]
  en = any(po.k .== k + 1) ? po.index[po.k .== k + 1][1] - 1 : st

  C = zeros(UInt8, m, k)
  v = zeros(Float64, k)
  n1s = zeros(Float64, m, k)

  M = -Inf

  p = zero(Float64)
  logprR = zero(Float64)
  logp = zeros(Float64, m)
  logpost = zero(Float64)

  for l in st:en
    R = po.partition[:, l]
    logprR = logdPriorRow(R, priorR)

    for g in 1:k
      v[g] = sum(R .== g)
      n1s[: , g] = sum(float(data[:, R .== g]), 2)
    end

    for c1 in cs, c2 in cs, c3 in cs, c4 in cs
      C[1, :] = c1
      C[2, :] = c2
      C[3, :] = c3
      C[4, :] = c4

      idx = C[:, 1] + 4 * collect(0:(m - 1))
      logp = priorC.logγ[C[:, 1]] + priorC.logω[C[:, 1]] +
             logmarglik(n1s[:, 1], v[1], priorC.A[idx], priorC.B[idx])

      for g in 2:k
        idx = C[:, g] + 4 * collect(0:(m - 1))
        logp += (priorC.logω[C[:, g]] +
                 logmarglik(n1s[:, g], v[g], priorC.A[idx], priorC.B[idx]))
      end

      logpost = logprR + sum(logp)

      if logpost > M
        M = logpost
      end

      p += exp(logpost - lumpp)
    end
  end

  (M, p)
end

function lognormconst(cs,
                      data::Array{UInt8, 2},
                      po::TestPartition,
                      γ::Array{Float64, 1},
                      r::Float64,
                      priorR::PriorRowPartition)
  # log unnormalized maximum posterior probability
  lumpp = -Inf

  println("Computing 'lumpp'...")
  for k in 1:size(data, 2)
    println("k = ", k)
    t1, t2 = computelognormconst(cs[k], k, 0.0, data, po, γ, r, priorR)

    if t1 > lumpp
      lumpp = t1
    end
  end
  println("Done.")

  # now that we know the maximum value, we can compute the logarithm of the
  # normalization constant
  p = zero(Float64)

  println("Computing 'p'...")
  for k in 1:size(data, 2)
    println("k = ", k)
    t1, t2 = computelognormconst(cs[k], k, lumpp, data, po, γ, r, priorR)
    p += t2
  end
  println("Done.")

  (log(p), lumpp)
end

function computeProbs(cs,
                      lp::Float64,
                      lumpp::Float64,
                      data::Array{UInt8, 2},
                      po::TestPartition,
                      γ::Array{Float64, 1},
                      r::Float64,
                      priorR::PriorRowPartition)
  m, n = size(data)

  p = zeros(Float64, div(n * (n - 1), 2))
  s = zeros(Float64, m , 3)

  u = falses(div(n * (n - 1), 2))

  logprR = zero(Float64)

  println("Computing probabilities...")
  for k in 1:(n - 1)
    println("k = ", k)
    priorC = AminoAcidPriorCol(data, k, γ, r)

    st = po.index[po.k .== k][1]
    en = po.index[po.k .== k + 1][1] - 1

    C = zeros(UInt8, m, k)
    v = zeros(Float64, k)
    n1s = zeros(Float64, m, k)

    for l in st:en
      R = po.partition[:, l]
      logprR = logdPriorRow(R, priorR)

      for g in 1:k
        v[g] = sum(R .== g)
        n1s[: , g] = sum(float(data[:, R .== g]), 2)
      end

      idx = zero(Int)
      for i in 1:(n - 1)
        for j in (i + 1):n
          u[idx += 1] = R[i] == R[j]
        end
      end

      for c1 in cs[k], c2 in cs[k], c3 in cs[k], c4 in cs[k]
        C[1, :] = c1
        C[2, :] = c2
        C[3, :] = c3
        C[4, :] = c4

        idx = C[:, 1] + 4 * collect(0:(m - 1))
        logp = priorC.logγ[C[:, 1]] + priorC.logω[C[:, 1]] +
               logmarglik(n1s[:, 1], v[1], priorC.A[idx], priorC.B[idx])

        for g in 2:k
          idx = C[:, g] + 4 * collect(0:(m - 1))
          logp += (priorC.logω[C[:, g]] +
                   logmarglik(n1s[:, g], v[g], priorC.A[idx], priorC.B[idx]))
        end

        tmp = exp(logprR + sum(logp) - lumpp)

        p[u] += tmp

        for b in 1:m
          if C[b, 1] == 0x01
            s[b, 1] += tmp
          elseif C[b, 1] == 0x02
            s[b, 2] += tmp
          else
            s[b, 3] += tmp
          end
        end
      end
    end
  end

  # no units are in the same cluster
  k = n
  println("k = ", k)

  R = collect(1:k)
  logprR = logdPriorRow(R, priorR)

  C = zeros(UInt8, m, k)
  priorC = AminoAcidPriorCol(data, k, γ, r)

  v = ones(Float64, k)
  n1s = float(data)

  for c1 in cs[k], c2 in cs[k], c3 in cs[k], c4 in cs[k]
    C[1, :] = c1
    C[2, :] = c2
    C[3, :] = c3
    C[4, :] = c4

    idx = C[:, 1] + 4 * collect(0:(m - 1))
    logp = priorC.logγ[C[:, 1]] + priorC.logω[C[:, 1]] +
           logmarglik(n1s[:, 1], v[1], priorC.A[idx], priorC.B[idx])

    for g in 2:k
      idx = C[:, g] + 4 * collect(0:(m - 1))
      logp += (priorC.logω[C[:, g]] +
               logmarglik(n1s[:, g], v[g], priorC.A[idx], priorC.B[idx]))
    end

    tmp = exp(logprR + sum(logp) - lumpp)

    for b in 1:m
      if C[b, 1] == 0x01
        s[b, 1] += tmp
      elseif C[b, 1] == 0x02
        s[b, 2] += tmp
      else
        s[b, 3] += tmp
      end
    end
  end
  println("Done.")

  (exp(log(p) - lp), exp(log(s) - lp))
end

x = AminoAcidData(fastafile)

settings = KSettings(outfile, T, burnin, tstep, op, α, θ, γ, r, λs1, λs2,
                     parawm, maxclust, maxunit, verbose, verbosestep)

priorR = EwensPitman(settings.α, settings.θ)
po = TestPartition(size(x.data, 2))
cs = ((UInt8[1], UInt8[2], UInt8[3], UInt8[4]),
      (UInt8[1; 1], UInt8[2; 2], UInt8[3; 3], UInt8[3; 4], UInt8[4; 3],
       UInt8[4; 4]),
      (UInt8[1; 1; 1], UInt8[2; 2; 2], UInt8[3; 3; 3], UInt8[3; 3; 4],
       UInt8[3; 4; 3], UInt8[4; 3; 3], UInt8[3; 4; 4], UInt8[4; 3; 4],
       UInt8[4; 4; 3], UInt8[4; 4; 4]),
      (UInt8[1; 1; 1; 1], UInt8[2; 2; 2; 2], UInt8[3; 3; 3; 3],
       UInt8[3; 3; 3; 4], UInt8[3; 3; 4; 3], UInt8[3; 4; 3; 3],
       UInt8[4; 3; 3; 3], UInt8[3; 3; 4; 4], UInt8[3; 4; 3; 4],
       UInt8[3; 4; 4; 3], UInt8[4; 3; 3; 4], UInt8[4; 3; 4; 3],
       UInt8[4; 4; 3; 3], UInt8[3; 4; 4; 4], UInt8[4; 3; 4; 4],
       UInt8[4; 4; 3; 4], UInt8[4; 4; 4; 3], UInt8[4; 4; 4; 4]),
      (UInt8[1; 1; 1; 1; 1], UInt8[2; 2; 2; 2; 2], UInt8[3; 3; 3; 3; 3],
       UInt8[3; 3; 3; 3; 4], UInt8[3; 3; 3; 4; 3], UInt8[3; 3; 4; 3; 3],
       UInt8[3; 4; 3; 3; 3], UInt8[4; 3; 3; 3; 3], UInt8[3; 3; 3; 4; 4],
       UInt8[3; 3; 4; 3; 4], UInt8[3; 3; 4; 4; 3], UInt8[3; 4; 3; 3; 4],
       UInt8[3; 4; 3; 4; 3], UInt8[3; 4; 4; 3; 3], UInt8[4; 3; 3; 3; 4],
       UInt8[4; 3; 3; 4; 3], UInt8[4; 3; 4; 3; 3], UInt8[4; 4; 3; 3; 3],
       UInt8[3; 3; 4; 4; 4], UInt8[3; 4; 3; 4; 4], UInt8[3; 4; 4; 3; 4],
       UInt8[3; 4; 4; 4; 3], UInt8[4; 3; 3; 4; 4], UInt8[4; 3; 4; 3; 4],
       UInt8[4; 4; 3; 3; 4], UInt8[4; 3; 4; 4; 3], UInt8[4; 4; 3; 4; 3],
       UInt8[4; 4; 4; 3; 3], UInt8[3; 4; 4; 4; 4], UInt8[4; 3; 4; 4; 4],
       UInt8[4; 4; 3; 4; 4], UInt8[4; 4; 4; 3; 4], UInt8[4; 4; 4; 4; 3],
       UInt8[4; 4; 4; 4; 4]),
      (UInt8[1; 1; 1; 1; 1; 1], UInt8[2; 2; 2; 2; 2; 2],
       UInt8[3; 3; 3; 3; 3; 3], UInt8[3; 3; 3; 3; 3; 4],
       UInt8[3; 3; 3; 3; 4; 3], UInt8[3; 3; 3; 4; 3; 3],
       UInt8[3; 3; 4; 3; 3; 3], UInt8[3; 4; 3; 3; 3; 3],
       UInt8[4; 3; 3; 3; 3; 3], UInt8[3; 3; 3; 3; 4; 4],
       UInt8[3; 3; 3; 4; 3; 4], UInt8[3; 3; 3; 4; 4; 3],
       UInt8[3; 3; 4; 3; 3; 4], UInt8[3; 3; 4; 3; 4; 3],
       UInt8[3; 3; 4; 4; 3; 3], UInt8[3; 4; 3; 3; 3; 4],
       UInt8[3; 4; 3; 3; 4; 3], UInt8[3; 4; 3; 4; 3; 3],
       UInt8[3; 4; 4; 3; 3; 3], UInt8[4; 3; 3; 3; 3; 4],
       UInt8[4; 3; 3; 3; 4; 3], UInt8[4; 3; 3; 4; 3; 3],
       UInt8[4; 3; 4; 3; 3; 3], UInt8[4; 4; 3; 3; 3; 3],
       UInt8[3; 3; 3; 4; 4; 4], UInt8[3; 3; 4; 3; 4; 4],
       UInt8[3; 3; 4; 4; 3; 4], UInt8[3; 3; 4; 4; 4; 3],
       UInt8[3; 4; 3; 3; 4; 4], UInt8[3; 4; 3; 4; 3; 4],
       UInt8[3; 4; 4; 3; 3; 4], UInt8[3; 4; 3; 4; 4; 3],
       UInt8[3; 4; 4; 3; 4; 3], UInt8[3; 4; 4; 4; 3; 3],
       UInt8[4; 3; 3; 3; 4; 4], UInt8[4; 3; 3; 4; 3; 4],
       UInt8[4; 3; 4; 3; 3; 4], UInt8[4; 4; 3; 3; 3; 4],
       UInt8[4; 3; 3; 4; 4; 3], UInt8[4; 3; 4; 3; 4; 3],
       UInt8[4; 4; 3; 3; 4; 3], UInt8[4; 3; 4; 4; 3; 3],
       UInt8[4; 4; 3; 4; 3; 3], UInt8[4; 4; 4; 3; 3; 3],
       UInt8[3; 3; 4; 4; 4; 4], UInt8[3; 4; 3; 4; 4; 4],
       UInt8[3; 4; 4; 3; 4; 4], UInt8[3; 4; 4; 4; 3; 4],
       UInt8[3; 4; 4; 4; 4; 3], UInt8[4; 3; 3; 4; 4; 4],
       UInt8[4; 3; 4; 3; 4; 4], UInt8[4; 3; 4; 4; 3; 4],
       UInt8[4; 3; 4; 4; 4; 3], UInt8[4; 4; 3; 3; 4; 4],
       UInt8[4; 4; 3; 4; 3; 4], UInt8[4; 4; 3; 4; 4; 3],
       UInt8[4; 4; 4; 3; 3; 4], UInt8[4; 4; 4; 3; 4; 3],
       UInt8[4; 4; 4; 4; 3; 3], UInt8[3; 4; 4; 4; 4; 4],
       UInt8[4; 3; 4; 4; 4; 4], UInt8[4; 4; 3; 4; 4; 4],
       UInt8[4; 4; 4; 3; 4; 4], UInt8[4; 4; 4; 4; 3; 4],
       UInt8[4; 4; 4; 4; 4; 3], UInt8[4; 4; 4; 4; 4; 4]))

# lp = 4.09068229126453619670655825757421553134918212890625
# lumpp = -21.0913924952027400649967603385448455810546875
# lc = lumpp + lp = -17.00071020393820475646862178109586238861083984375
# lmpp = lumpp - lc = lumpp - lumpp - lp = - lp
# lc + lmpp = lumpp + lp - lp = lumpp
lp, lumpp = lognormconst(cs, x.data, po, settings.γ, settings.r, priorR)
p, s = computeProbs(cs, lp, lumpp, x.data, po, settings.γ, settings.r, priorR)
=#

p = [0.5647226958512603367523752240231260657310485839843750000;
     0.3859601712857648747601047034549992531538009643554687500;
     0.3859601712857648747601047034549992531538009643554687500;
     0.3877266490089570361021742428420111536979675292968750000;
     0.2495816768360757109679326504192431457340717315673828125;
     0.3859601712857648747601047034549992531538009643554687500;
     0.3859601712857648747601047034549992531538009643554687500;
     0.3877266490089570361021742428420111536979675292968750000;
     0.2495816768360757109679326504192431457340717315673828125;
     0.5647226958512601147077702989918179810047149658203125000;
     0.2495816768360775983470745131853618659079074859619140625;
     0.3877266490089560369014520802011247724294662475585937500;
     0.2495816768360764881240498880288214422762393951416015625;
     0.3877266490089572026356279366154922172427177429199218750;
     0.3775410735825202035442771375528536736965179443359375000]

s = reshape([0.646616210995535456440563848445890471339225769042968750000;
             0.646616210995536011552076161024160683155059814453125000000;
             0.685422396679224998905510801705531775951385498046875000000;
             0.681017613877477279160643774957861751317977905273437500000;
             0.321786066991730956843298372405115514993667602539062500000;
             0.321786066991730346220634828569018281996250152587890625000;
             0.288655590615555790456880913552595302462577819824218750000;
             0.286366638965632747115819256578106433153152465820312500000;
             0.031597722012356437015778709564983728341758251190185546875;
             0.031597722012359129306613425569594255648553371429443359375;
             0.025922012704874257404963344697534921579062938690185546875;
             0.032615747156614409429931100703470292501151561737060546875],
            (4, 3))

kpax3aa(AminoAcidData(fastafile), partition, outfile, T, burnin=burnin,
        tstep=tstep, op=op, α=α, θ=θ, γ=γ, r=r, λs1=λs1, λs2=λs2, parawm=parawm,
        maxclust=maxclust, maxunit=maxunit, verbose=verbose,
        verbosestep=verbosestep)
