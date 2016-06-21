# This file is part of Kpax3. License is MIT.

# TODO: test normalizepartition, modifypartition, modifymerge, modifysplit,
#       modifyscramble

function test_partition_functions_initializepartition()
  ifile = "data/proper_aa.fasta"
  ofile = "../build/test"

  settings = KSettings(ifile, ofile)

  R = initializepartition(settings)

  @test isa(R, Vector{Int})
  @test minimum(R) >= 1
  @test maximum(R) <= 6

  R = initializepartition(settings, kset=1:6)

  @test isa(R, Vector{Int})
  @test minimum(R) >= 1
  @test maximum(R) <= 6

  nothing
end

test_partition_functions_initializepartition()

function test_partition_functions_encodepartition()
  R = zeros(Int, 6)
  S = [3; 3; 3; 2; 2; 1]

  encodepartition!(R, S)

  @test R == [1; 2; 3; 22; 23; 36]

  nothing
end

test_partition_functions_encodepartition()

function test_partition_functions_decodepartition()
  R = [1; 2; 3; 22; 23; 36]

  decodepartition!(R)

  @test R == [1; 1; 1; 4; 4; 6]

  nothing
end

test_partition_functions_decodepartition()
