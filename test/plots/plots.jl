# This file is part of Kpax3. License is MIT.

function test_plotk()
  (k, pk) = Kpax3.readposteriork("../build/mcmc_6")
  p = Kpax3.plotk(k, pk, xticks=[1; 2; 3; 4; 5; 6], width=800, height=600)

  Plots.png(p, "../build/MCMC_posterior_k.png")

  nothing
end

test_plotk()

function test_plotC()
  (site, aa, freq, C) = Kpax3.readposteriorC("../build/mcmc_6")

  p = Kpax3.plotC(site, freq, C, width=800, height=600)

  Plots.png(p, "../build/MCMC_posterior_column_classifier.png")

  nothing
end

test_plotC()

function test_plotD()
  settings = Kpax3.KSettings("data/mcmc_6.fasta", "../build/tmp")
  x = Kpax3.AminoAcidData(settings)
  state = Kpax3.optimumstate(x, "data/mcmc_6.csv", settings)

  p = Kpax3.plotD(x, state, clusterorder=[4; 2; 1; 3],
                  clusterlabel=["d"; "b"; "a"; "c"], width=800, height=600)

  Plots.png(p, "../build/MCMC_posterior_dataset.png")

  nothing
end

test_plotD()

function test_plotP()
  (id, P) = Kpax3.readposteriorP("../build/mcmc_6")
  R = Kpax3.normalizepartition("data/mcmc_6.csv", id)

  p = Kpax3.plotP(R, P, clusterorder=[4; 2; 1; 3],
                  clusterlabel=["d"; "b"; "a"; "c"], width=800.0, height=800.0)

  Plots.png(p, "../build/MCMC_posterior_probability_matrix.png")

  nothing
end

test_plotP()

function test_plottrace()
  (entropy_R, avgd_R) = Kpax3.traceR("../build/mcmc_6", maxlag=50)
  (entropy_C, avgd_C) = Kpax3.traceC("../build/mcmc_6", maxlag=50)

  p = Kpax3.plottrace(entropy_R, maxlag=50, main="Trace",
                      width=800, height=600)
  q = Kpax3.plottrace(entropy_C, maxlag=50, main="Trace",
                      width=800, height=600)

  Plots.png(p, "../build/MCMC_posterior_trace_R.png")
  Plots.png(q, "../build/MCMC_posterior_trace_C.png")

  nothing
end

test_plottrace()

function test_plotdensity()
  (entropy_R, avgd_R) = Kpax3.traceR("../build/mcmc_6", maxlag=50)
  (entropy_C, avgd_C) = Kpax3.traceC("../build/mcmc_6", maxlag=50)

  p = Kpax3.plotdensity(entropy_R, maxlag=50, main="Density",
                        width=800, height=600)
  q = Kpax3.plotdensity(entropy_C, maxlag=50, main="Density",
                        width=800, height=600)

  Plots.png(p, "../build/MCMC_posterior_density_R.png")
  Plots.png(q, "../build/MCMC_posterior_density_C.png")

  nothing
end

test_plotdensity()

function test_plotjump()
  (entropy_R, avgd_R) = Kpax3.traceR("../build/mcmc_6", maxlag=50)
  (entropy_C, avgd_C) = Kpax3.traceC("../build/mcmc_6", maxlag=50)

  p = Kpax3.plotjump(avgd_R, main="Jump distance", width=800, height=600)
  q = Kpax3.plotjump(avgd_C, main="Jump distance", width=800, height=600)

  Plots.png(p, "../build/MCMC_posterior_jump_R.png")
  Plots.png(q, "../build/MCMC_posterior_jump_C.png")

  nothing
end

test_plotjump()

function test_plotdgn()
  (entropy_R, avgd_R) = Kpax3.traceR("../build/mcmc_6", maxlag=50)
  (entropy_C, avgd_C) = Kpax3.traceC("../build/mcmc_6", maxlag=50)

  p = Kpax3.plotdgn(entropy_R, avgd_R, main="Sample Partition", maxlag=50,
                    width=800, height=600)
  q = Kpax3.plotdgn(entropy_C, avgd_C, main="Sample Partition", maxlag=50,
                    width=800, height=600)

  Plots.png(p, "../build/MCMC_posterior_diagnosis_R.png")
  Plots.png(q, "../build/MCMC_posterior_diagnosis_C.png")

  nothing
end

test_plotdgn()
