# This file is part of Kpax3. License is MIT.

v = 4
unit = [1; 2; 3; 4]
n1s = Float64[0; 4; 1; 3; 2; 1; 3; 1; 1; 2; 1; 0; 3; 1; 2; 2; 3; 1]

cluster = KCluster(v, unit, n1s)

# TODO: proper testing
