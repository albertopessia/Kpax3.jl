# This file is part of Kpax3. License is MIT.

m = 18
n = 6
maxclust = 1
maxunit = 1

support = KSupport(m, n, maxclust, maxunit)

@test support.m == m
@test support.n == n

@test support.vi == 0
@test support.ni == zeros(Float64, m)
@test support.ui == zeros(Int, maxunit)
@test isa(support.wi, KWeight)
@test support.wi.c == zeros(Float64, m)
@test support.wi.w == zeros(Float64, 4, m)
@test support.wi.z == zeros(Float64, 4, m)

@test support.vj == 0
@test support.nj == zeros(Float64, m)
@test support.uj == zeros(Int, maxunit)
@test isa(support.wj, KWeight)
@test support.wj.c == zeros(Float64, m)
@test support.wj.w == zeros(Float64, 4, m)
@test support.wj.z == zeros(Float64, 4, m)

@test support.tmp == zeros(Float64, 4)

@test support.cl == zeros(Int, n)
@test support.k == 0

@test isa(support.oi, KOffspring)
@test support.oi.R == zeros(Int, n)
@test support.oi.v == zeros(Int, n)

@test isa(support.oj, KOffspring)
@test support.oj.R == zeros(Int, n)
@test support.oj.v == zeros(Int, n)

@test support.C == zeros(UInt8, maxclust, m)

@test support.logpC == zeros(Float64, 2)
@test support.lograR == 0.0
@test support.loglik == 0.0
