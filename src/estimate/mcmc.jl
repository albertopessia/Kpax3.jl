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

  lossold = Inf
  lossnew = Inf
  tmp = 0
  for k in 1:n
    if pk[k] > 0.0
      estimate = kmedoids(D, k).assignments
      tmp = 1
      while tmp <= 100
        lossnew = loss_binder(estimate, pR)

        if lossnew < lossold
          lossold = lossnew
          copy!(R, estimate)
        end

        estimate = kmedoids(D, k).assignments
        tmp += 1
      end
    end
  end

  R
end

function kpax3estimate(x::AminoAcidData,
                       settings::KSettings)
  R = kpax3Restimate(settings.ofile)

  (m, n) = size(x.data)

  fp = open(settings.ofile, "r")
  tmp = zeros(Float64, 6)
  read!(fp, tmp)
  close(fp)

  α = tmp[1]
  θ = tmp[2]
  γ = [tmp[3]; tmp[4]; tmp[5]]
  r = tmp[6]

  op = copy(values(settings.op))
  (λs1, λs2) = params(settings.distws)

  settings = KSettings(settings.ifile, settings.ofile, α, θ, γ, r,
                       settings.maxclust, settings.maxunit, settings.verbose,
                       settings.verbosestep, settings.popsize, settings.xrate,
                       settings.mrate, settings.T, settings.burnin,
                       settings.tstep, WeightVec(op), Beta(λs1, λs2),
                       settings.parawm)

  k = maximum(R)

  priorR = EwensPitman(α, θ)
  priorC = AminoAcidPriorCol(x.data, γ, r)

  AminoAcidState(x.data, R, priorR, priorC, settings)
end
