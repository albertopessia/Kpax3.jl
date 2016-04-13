# This file is part of Kpax3. License is MIT.

function kpax3Restimate(infile::AbstractString)
  (pk, pR, pS) = readresults(infile)

  n = length(pk)

  D = zeros(Float64, n, n)
  idx = 0
  for j in 1:(n - 1), i in (j + 1):n
    idx += 1
    D[i, j] = D[j, i] = 1.0 - pR[idx]
  end

  R = zeros(Int, n)

  lossold = Inf
  lossnew = Inf
  tmp = 0
  for k in 1:n
    if pk[k] > 0.0
      estimate = kmedoids(D, k; maxiter=1000)
      tmp = 1
      while tmp <= 100
        lossnew = loss_binder(estimate.assignments, pR)

        if lossnew < lossold
          lossold = lossnew
          copy!(R, estimate.assignments)
        end

        estimate = kmedoids(D, k; maxiter=1000)
        tmp += 1
      end
    end
  end

  R
end

function kpax3estimate(x::AminoAcidData,
                       infile::AbstractString,
                       settings::KSettings)
  R = kpax3Restimate(infile)

  (m, n) = size(x.data)

  fp = open(infile, "r")
  tmp = zeros(Float64, 6)
  read!(fp, tmp)
  close(fp)

  α = tmp[1]
  θ = tmp[2]
  γ = [tmp[3]; tmp[4]; tmp[5]]
  r = tmp[6]

  settings = KSettings(settings.fpath, settings.T, settings.burnin,
                       settings.tstep, settings.op, α, θ, γ, r, settings.distws,
                       settings.parawm, settings.maxclust, settings.maxunit,
                       settings.verbose, settings.verbosestep)

  k = maximum(R)

  priorR = EwensPitman(α, θ)
  priorC = AminoAcidPriorCol(x.data, k, γ, r)

  AminoAcidState(x.data, R, priorR, priorC, settings)
end
