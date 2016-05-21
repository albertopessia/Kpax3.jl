# This file is part of Kpax3. License is MIT.

ifile = "data/proper_aa.fasta"
ofile = "../build/test.bin"

settings = KSettings(ifile, ofile)

R = initializepartition(settings)

@test isa(R, Vector{Int})
@test minimum(R) >= 1
@test maximum(R) <= 6

R = initializepartition(settings, kset=1:6)

@test isa(R, Vector{Int})
@test minimum(R) >= 1
@test maximum(R) <= 6

# TODO: test normalizepartition functions

R = zeros(Int, 6)
S = [3; 3; 3; 2; 2; 1]

encodepartition!(R, S)

@test R == [1; 2; 3; 22; 23; 36]
