# This file is part of Kpax3. License is MIT.

# TODO: test normalizepartition functions

R = zeros(Int, 6)
S = [3; 3; 3; 2; 2; 1]

encodepartition!(R, S)

@test R == [1; 2; 3; 22; 23; 36]
