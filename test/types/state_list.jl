# This file is part of Kpax3. License is MIT.

ifile = "data/proper_aa.fasta"
ofile = "../build/test.bin"

settings = KSettings(ifile, ofile)

x = AminoAcidData(settings)

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(x.data, settings.γ, settings.r)

R = zeros(Int, size(x.data, 2))

slist = AminoAcidStateList(x.data, [3; 3; 3; 1; 1; 2], priorR, priorC, settings)

@test isa(slist.state, Vector{AminoAcidState})
@test isa(slist.logpp, Vector{Float64})
@test isa(slist.rank, Vector{Int})
@test all(slist.logpp .< 0.0)
@test length(unique(slist.logpp)) > 1
@test all(0 .< slist.rank .< settings.popsize + 1)
@test length(unique(slist.rank)) == settings.popsize
@test all(diff(slist.logpp[slist.rank]) .<= 0)

@test slist.state[1].R == [1; 2; 3; 22; 23; 36]
#=
for i in 1:settings.popsize
  t = AminoAcidState(x.data, slist.state[i].R, priorR, priorC, settings)
  tlp = t.logpR + t.logpC[1] + t.loglik

  l = t.cl[1:t.k]

  @test slist.state[i].R == encodepartition(t.R)
  @test slist.state[i].C == t.C
  @test slist.state[i].emptycluster == t.emptycluster
  @test slist.state[i].cl == t.cl
  @test slist.state[i].k == t.k
  @test slist.state[i].v == t.v
  @test slist.state[i].n1s == t.n1s
  for g in l
    @test slist.state[i].unit[g][1:slist.state[i].v[g]] == t.unit[g][1:t.v[g]]
  end
  @test slist.state[i].logpR == t.logpR
  @test slist.state[i].logpC == t.logpC
  @test slist.state[i].loglik == t.loglik

  @test_approx_eq_eps slist.logpp[i] tlp ε
end
=#
(m, n) = size(x.data)

D = zeros(Float64, n, n)
for j in 1:(n - 1), i in (j + 1):n
  D[i, j] = D[j, i] = sum(x.data[:, j] .!= x.data[:, i]) / m
end

slist = AminoAcidStateList(x.data, D, 1:n, priorR, priorC, settings)

@test isa(slist.state, Vector{AminoAcidState})
@test isa(slist.logpp, Vector{Float64})
@test isa(slist.rank, Vector{Int})
@test all(slist.logpp .< 0.0)
@test length(unique(slist.logpp)) > 1
@test all(0 .< slist.rank .< settings.popsize + 1)
@test length(unique(slist.rank)) == settings.popsize
@test all(diff(slist.logpp[slist.rank]) .<= 0)
#=
for i in 1:settings.popsize
  t = AminoAcidState(x.data, slist.state[i].R, priorR, priorC, settings)
  tlp = t.logpR + t.logpC[1] + t.loglik

  l = t.cl[1:t.k]

  @test slist.state[i].R == t.R
  @test slist.state[i].C == t.C
  @test slist.state[i].emptycluster == t.emptycluster
  @test slist.state[i].cl == t.cl
  @test slist.state[i].k == t.k
  @test slist.state[i].v == t.v
  @test slist.state[i].n1s == t.n1s
  for g in l
    @test slist.state[i].unit[g][1:slist.state[i].v[g]] == t.unit[g][1:t.v[g]]
  end
  @test slist.state[i].logpR == t.logpR
  @test slist.state[i].logpC == t.logpC
  @test slist.state[i].loglik == t.loglik

  @test_approx_eq_eps slist.logpp[i] tlp ε
end
=#
