# This file is part of Kpax3. License is MIT.

function kpax3Restimate(ifile::AbstractString)
  (pk, pR, pS) = processchain(ifile)

  n = length(pk)

  D = zeros(Float64, n, n)
  idx = 0
  for j in 1:(n - 1), i in (j + 1):n
    idx += 1
    D[i, j] = D[j, i] = 1.0 - pR[idx]
  end

  R = zeros(Int, n)
  estimate = ones(Int, n)

  lossold = Inf
  lossnew = Inf
  niter = 0
  for k in 1:n
    if pk[k] > 0.0
      try
        copy!(estimate, kmedoids(D, k).assignments)
      catch
        sample!(1:k, estimate, replace=true)
        estimate[sample(1:n, k, replace=false)] = collect(1:k)
      end

      niter = 0
      while niter < 100
        lossnew = loss_binder(estimate, pR)

        if lossnew < lossold
          lossold = lossnew
          copy!(R, estimate)
        end

        try
          copy!(estimate, kmedoids(D, k).assignments)
        catch
          sample!(1:k, estimate, replace=true)
          estimate[sample(1:n, k, replace=false)] = collect(1:k)
        end

        niter += 1
      end
    end
  end

  R
end

function kpax3estimate(x::AminoAcidData,
                       settings::KSettings)
  R = kpax3Restimate(settings.ofile)

  fpS = open(string(settings.ofile, "_settings.bin"), "r")

  tmp1 = zeros(Int, 3)
  read!(fpS, tmp1)

  n = tmp1[1]
  m = tmp1[2]
  N = tmp1[3]

  tmp2 = zeros(Float64, 6)
  read!(fpS, tmp2)

  α = tmp2[1]
  θ = tmp2[2]
  γ = [tmp2[3]; tmp2[4]; tmp2[5]]
  r = tmp2[6]

  close(fpS)

  op = copy(values(settings.op))

  settings = KSettings(settings.ifile, settings.ofile, settings.protein,
                       settings.miss, settings.l, α, θ, γ, r, settings.maxclust,
                       settings.maxunit, settings.verbose, settings.verbosestep,
                       settings.popsize, settings.maxiter, settings.maxgap,
                       settings.xrate, settings.mrate, settings.T,
                       settings.burnin, settings.tstep, WeightVec(op))

  k = maximum(R)

  priorR = EwensPitman(α, θ)
  priorC = AminoAcidPriorCol(x.data, γ, r)

  AminoAcidState(x.data, R, priorR, priorC, settings)
end
