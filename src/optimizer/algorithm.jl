# This file is part of Kpax3. License is MIT.

function kpax3ga!(x::AminoAcidData,
                  population::AminoAcidStateList,
                  priorR::PriorRowPartition,
                  priorC::PriorColPartition,
                  settings::KSettings,
                  support::KSupport)
  if settings.verbose
    @printf("Initializing Genetic Algorithm... ")
  end

  # check if we can write to the backup file
  fp = open(settings.ofile, "w")
  close(fp)

  # elitism
  Nelite = max(1, ceil(Int, settings.popsize * 0.2))

  # randomization
  Nrand = Nelite + floor(Int, settings.popsize * 0.1)

  if ((settings.popsize - Nrand) % 2) != 0
    Nrand += 1
  end

  # initialize variables
  R = zeros(Int, support.n)

  i = 1
  maxk = 0
  idx = 0
  for i in 1:settings.popsize
    if population.state[i].k > maxk
      maxk = population.state[i].k
      idx = i
    end
  end

  beststate = copystate(population.state[population.rank[1]])

  newpopulation = AminoAcidStateList(settings.popsize, population.state[idx])

  if settings.verbose
    @printf("done.\nStochastic optimization will now begin.\n")
  end

  iter = 0
  gap = 0
  keepgoing = true
  while keepgoing
    # copy the first Nelite best solutions without changing them
    i = 1
    while i <= Nelite
      copystate!(newpopulation.state[i], population.state[population.rank[i]])
      newpopulation.logpp[i] = population.logpp[population.rank[i]]
      i += 1
    end

    # now create Nrand random solutions starting from the best one
    while i <= Nrand
      copy!(R, newpopulation.state[1].R)
      modifypartition!(R, newpopulation.state[1].k)

      newpopulation.state[i] = AminoAcidState(x.data, R, priorR, priorC,
                                              settings)
      newpopulation.logpp[i] = newpopulation.state[i].logpp

      i += 1
    end

    while i <= settings.popsize
      (i1, i2) = selection(population.logpp)

      if rand() <= settings.xrate
        crossover!(population.state[i1].R, population.state[i2].R, support)
      else
        copy!(support.oi.R, population.state[i1].R)
        copy!(support.oi.v, population.state[i1].v)

        copy!(support.oj.R, population.state[i2].R)
        copy!(support.oj.v, population.state[i2].v)
      end

      mutation!(support.oi, settings.mrate)
      mutation!(support.oj, settings.mrate)

      updatestate!(newpopulation.state[i], x.data, support.oi.R, priorR, priorC,
                   settings)
      updatestate!(newpopulation.state[i + 1], x.data, support.oj.R, priorR,
                   priorC, settings)

      newpopulation.logpp[i] = newpopulation.state[i].logpp
      newpopulation.logpp[i + 1] = newpopulation.state[i + 1].logpp

      i += 2
    end

    sortperm!(newpopulation.rank, newpopulation.logpp, rev=true,
              initialized=true)

    copystatelist!(population, newpopulation, settings.popsize)

    if population.logpp[population.rank[1]] >  beststate.logpp
      copystate!(beststate, population.state[population.rank[1]])
      gap = 0

      if settings.verbose
        @printf("Found a better solution! ")
        @printf("Log-posterior (plus a constant) for %d clusters: %.4f.\n",
                beststate.k, beststate.logpp)
      end
    else
      gap += 1
      keepgoing = (gap <= settings.maxgap)
    end

    iter += 1

    if iter < settings.maxiter
      if (iter % settings.verbosestep == 0)
        fp = open(settings.ofile, "w")
        for i in 1:support.n
          write(fp, "\"$(x.id[i])\",$(beststate.R[i])\n")
        end
        close(fp)

        if settings.verbose
          println("Step ", iter, " done.")
        end
      end
    else
      keepgoing = false
    end
  end

  beststate
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
  (m, n) = size(x.data)

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r,
                             maxclust=max(kset[end], settings.maxclust))

  population = AminoAcidStateList(x.data, D, kset, priorR, priorC, settings)

  maxunit = min(n, max(maximum(population.state[population.rank[1]].v),
                       settings.maxunit))

  support = KSupport(m, n, kset[end], maxunit)

  kpax3ga!(x, population, priorR, priorC, settings, support)
end

function kpax3ga(x::AminoAcidData,
                 partition,
                 settings::KSettings)
  (m, n) = size(x.data)

  R = normalizepartition(partition, x.id)
  k = maximum(R)

  maxclust = min(n, max(k, settings.maxclust))

  priorR = EwensPitman(settings.α, settings.θ)
  priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r, maxclust=maxclust)

  population = AminoAcidStateList(x.data, R, priorR, priorC, settings)

  maxunit = min(n, max(maximum(population.state[population.rank[1]].v),
                       settings.maxunit))

  support = KSupport(m, n, maxclust, maxunit)

  kpax3ga!(x, population, priorR, priorC, settings, support)
end
