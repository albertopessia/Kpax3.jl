# This file is part of Kpax3. License is MIT.

function kpax3mcmc(settings::KSettings)
  R = normalizepartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(k, settings.maxclust))

  support = KSupport(size(x.data, 1), size(x.data, 2), settings.maxclust,
                     settings.maxunit)

  state = AminoAcidState(x.data, R, priorR, priorC, settings)

  kpax3mcmc!(x.data, priorR, priorC, settings, support, state)

  nothing
end

function kpax3mcmc(x::AminoAcidData,
                   partition,
                   settings::KSettings)
  R = normalizepartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(k, settings.maxclust))

  support = KSupport(size(x.data, 1), size(x.data, 2), settings.maxclust,
                     settings.maxunit)

  state = AminoAcidState(x.data, R, priorR, priorC, settings)

  kpax3mcmc!(x.data, priorR, priorC, settings, support, state)

  nothing
end

function kpax3ga(settings::KSettings;
                 kset::UnitRange{Int}=1:0)
  if settings.verbose
    @printf("Computing pairwise distances... ")
  end

  tmp = zeros(UInt8, length(settings.miss))
  idx = 0
  for c in 1:length(settings.miss)
    if settings.miss[c] != UInt8('-')
      idx += 1
      tmp[idx] = settings.miss[c]
    end
  end

  miss = if idx > 0
           copy!(zeros(UInt8, idx), 1, tmp, 1, idx)
         else
           zeros(UInt8, 1)
         end

  (data, id, ref) = readfasta(settings.ifile, settings.protein, miss,
                              settings.l, false, 0)

  n = size(data, 2)

  d = if settings.protein
        distaamtn84(data, ref)
      else
        distntmtn93(data, ref)
      end

  D = zeros(Float64, n, n)
  idx = 1
  for j in 1:(n - 1), i in (j + 1):n
    D[i, j] = D[j, i] = d[idx]
    idx += 1
  end

  if settings.verbose
    @printf("done.\n")
  end

  # expected number of cluster approximately between cbrt(n) and sqrt(n)
  g = ceil(Int, n^(2 / 5))

  kset = if length(kset) == 0
           max(1, g - 20):min(n, g + 20)
         elseif kset[1] > 0
           if kset[end] > n
             if kset[1] == 1
               2:n
             else
               kset[1]:n
             end
           elseif kset[1] == 1
             2:kset[end]
           end
         else
           throw(KDomainError("First element of 'kset' is less than one."))
         end

  x = AminoAcidData(settings)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(kset[end], settings.maxclust))

  population = AminoAcidStateList(x.data, D, kset, priorR, priorC, settings)

  kpax3ga!(x.data, population, priorR, priorC, settings)

  nothing
end

function kpax3ga(x::AminoAcidData,
                 partition,
                 settings::KSettings)
  R = normalizepartition(partition, x.id)
  k = maximum(R)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(k, settings.maxclust))

  population = AminoAcidStateList(x.data, R, priorR, priorC, settings)

  kpax3ga!(x.data, population, priorR, priorC, settings)

  nothing
end
