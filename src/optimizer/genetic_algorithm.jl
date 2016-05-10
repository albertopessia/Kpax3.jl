# This file is part of Kpax3. License is MIT.

function kpax3ga!(data::Matrix{UInt8},
                  population::AminoAcidStateList,
                  priorR::PriorRowPartition,
                  priorC::PriorColPartition,
                  settings::KSettings)
  # check if we can write to the backup file
  fp = open(settings.ofile, "w")
  close(fp)

  # elitism
  N = max(1, ceil(Int, settings.popsize * 0.2))

  # initialize new population
  newpopulation = AminoAcidStateList(settings.popsize, population.state[1])

  iter = 0
  gap = 0
  keepgoing = true
  while keepgoing
  end

  nothing
end
