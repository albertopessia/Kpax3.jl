# This file is part of Kpax3. License is MIT.

function test_write_results()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"
  partition = "data/proper_aa_partition.csv"

  settings = Kpax3.KSettings(ifile, ofile, gamma=[0.4; 0.35; 0.25])

  x = Kpax3.AminoAcidData(settings)

  R = Kpax3.normalizepartition(partition, x.id)
  k = maximum(R)

  priorR = Kpax3.EwensPitman(settings.α, settings.θ)
  priorC = Kpax3.AminoAcidPriorCol(x.data, settings.γ, settings.r)

  state = Kpax3.AminoAcidState(x.data, R, priorR, priorC, settings)

  Kpax3.writeresults(x, state, "../build/test_results", what=4)

  y1 = readstring("../build/test_results_partition.csv")
  y2 = parse(Float64, strip(readstring("../build/test_results_logposterior_value.txt")))
  y3 = readstring("../build/test_results_attributes.csv")
  y4 = readstring("../build/test_results_characteristic.csv")
  y5 = readstring("../build/test_results_dataset.txt")

  @test y1 == readstring(partition)
  @test_approx_eq_eps y2 state.logpp ε
  @test y3 == readstring("data/proper_aa_attributes.csv")
  @test y4 == readstring("data/proper_aa_characteristic.csv")
  @test y5 == readstring("data/proper_aa_dataset.txt")

  nothing
end

test_write_results()
