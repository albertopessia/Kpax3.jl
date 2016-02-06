# This file is part of Kpax3. License is MIT.

ε = 2.0e-14

data = UInt8[0x00 0x00 0x00 0x00 0x00 0x01;
             0x01 0x01 0x01 0x01 0x01 0x00;
             0x00 0x00 0x01 0x00 0x01 0x01;
             0x01 0x01 0x00 0x01 0x00 0x00;
             0x01 0x01 0x00 0x00 0x00 0x00;
             0x00 0x00 0x00 0x01 0x01 0x00;
             0x01 0x01 0x01 0x00 0x00 0x00;
             0x00 0x00 0x00 0x01 0x01 0x01;
             0x00 0x00 0x01 0x00 0x00 0x00;
             0x01 0x00 0x00 0x01 0x00 0x01;
             0x00 0x01 0x00 0x00 0x01 0x00;
             0x00 0x00 0x00 0x00 0x00 0x01;
             0x01 0x01 0x01 0x00 0x00 0x00;
             0x00 0x00 0x00 0x01 0x01 0x00;
             0x01 0x01 0x00 0x00 0x01 0x01;
             0x00 0x00 0x01 0x01 0x00 0x00;
             0x01 0x01 0x00 0x01 0x00 0x00;
             0x00 0x00 0x01 0x00 0x01 0x01]

m, n = size(data)

settings = KSettings("typesmcmc.bin", 1, 0, 1, [1.0; 0.0; 0.0], 0.0, 1.0,
                     [0.6; 0.35; 0.05], 135.0, 1.0, 1.0, 5.0, 100, 1, true, 1)

R = [13; 13; 13; 13; 42; 42]

k = length(unique(R))

filledcluster = falses(settings.maxclust)
filledcluster[unique(R)] = true

cl = find(filledcluster)

v = zeros(Int, settings.maxclust)
v[cl[1]] = 4
v[cl[2]] = 2

n1s = zeros(Float64, settings.maxclust, m)
n1s[cl[1], :] = [0.0; 4.0; 1.0; 3.0; 2.0; 1.0; 3.0; 1.0; 1.0; 2.0; 1.0; 0.0;
                 3.0; 1.0; 2.0; 2.0; 3.0; 1.0]
n1s[cl[2], :] = [1.0; 1.0; 2.0; 0.0; 0.0; 1.0; 0.0; 2.0; 0.0; 1.0; 1.0; 1.0;
                 0.0; 1.0; 2.0; 0.0; 0.0; 2.0]

α = 0.0
θ = 1.0
priorR = EwensPitman(α, θ)

g = [0.6; 0.35; 0.05]
r = log(0.001) / log(0.95)
priorC = AminoAcidPriorCol(data, k, g, r)

mcmcobj = AminoAcidMCMC(data, R, priorR, priorC, settings)

@test mcmcobj.R == R

@test isa(mcmcobj.C, Array{UInt8, 2})
@test size(mcmcobj.C, 1) == settings.maxclust
@test size(mcmcobj.C, 2) == m

@test mcmcobj.filledcluster == filledcluster
@test mcmcobj.cl == cl

@test mcmcobj.v == v
@test mcmcobj.n1s == n1s
@test mcmcobj.unit == Vector{Int}[sum(R .== g) > 0 ? find(R .== g) : [0]
                                  for g in 1:settings.maxclust]

@test mcmcobj.logpR == logdPriorRow(n, k, v, priorR)
@test_approx_eq_eps mcmcobj.logpC[1] logpriorC(mcmcobj.C, mcmcobj.cl,
                                               priorC.logγ, priorC.logω) eps()
@test_approx_eq_eps mcmcobj.logpC[2] logcondpostC(mcmcobj.C, mcmcobj.cl,
                                                  mcmcobj.v, mcmcobj.n1s,
                                                  priorC.logω, priorC) eps()

# is linearidx approach correct?
C = UInt8[1 1 1 4 2 1 4 3 1 1 1 1 1 1 1 1 1 2;
          1 1 1 3 2 1 3 4 1 1 1 1 1 1 1 1 1 2]

A = [priorC.A[1,  1] priorC.A[1,  1];
     priorC.A[1,  2] priorC.A[1,  2];
     priorC.A[1,  3] priorC.A[1,  3];
     priorC.A[4,  4] priorC.A[3,  4];
     priorC.A[2,  5] priorC.A[2,  5];
     priorC.A[1,  6] priorC.A[1,  6];
     priorC.A[4,  7] priorC.A[3,  7];
     priorC.A[3,  8] priorC.A[4,  8];
     priorC.A[1,  9] priorC.A[1,  9];
     priorC.A[1, 10] priorC.A[1, 10];
     priorC.A[1, 11] priorC.A[1, 11];
     priorC.A[1, 12] priorC.A[1, 12];
     priorC.A[1, 13] priorC.A[1, 13];
     priorC.A[1, 14] priorC.A[1, 14];
     priorC.A[1, 15] priorC.A[1, 15];
     priorC.A[1, 16] priorC.A[1, 16];
     priorC.A[1, 17] priorC.A[1, 17];
     priorC.A[2, 18] priorC.A[2, 18]]

B = [priorC.B[1,  1] priorC.B[1,  1];
     priorC.B[1,  2] priorC.B[1,  2];
     priorC.B[1,  3] priorC.B[1,  3];
     priorC.B[4,  4] priorC.B[3,  4];
     priorC.B[2,  5] priorC.B[2,  5];
     priorC.B[1,  6] priorC.B[1,  6];
     priorC.B[4,  7] priorC.B[3,  7];
     priorC.B[3,  8] priorC.B[4,  8];
     priorC.B[1,  9] priorC.B[1,  9];
     priorC.B[1, 10] priorC.B[1, 10];
     priorC.B[1, 11] priorC.B[1, 11];
     priorC.B[1, 12] priorC.B[1, 12];
     priorC.B[1, 13] priorC.B[1, 13];
     priorC.B[1, 14] priorC.B[1, 14];
     priorC.B[1, 15] priorC.B[1, 15];
     priorC.B[1, 16] priorC.B[1, 16];
     priorC.B[1, 17] priorC.B[1, 17];
     priorC.B[2, 18] priorC.B[2, 18]]

loglik = zeros(Float64, 2)

linearidx = [(C[1, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[1] = sum(logmarglik(vec(n1s[cl[1], :]), v[cl[1]], priorC.A[linearidx],
                           priorC.B[linearidx]))

linearidx = [(C[2, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[2] = sum(logmarglik(vec(n1s[cl[2], :]), v[cl[2]], priorC.A[linearidx],
                           priorC.B[linearidx]))

@test loglik[1] == sum(logmarglik(vec(n1s[cl[1], :]), v[cl[1]], A[:, 1],
                                  B[:, 1]))
@test loglik[2] == sum(logmarglik(vec(n1s[cl[2], :]), v[cl[2]], A[:, 2],
                                  B[:, 2]))

loglik = zeros(Float64, 2)
linearidx = [(mcmcobj.C[cl[1], b] + 4 * (b - 1))::Int for b in 1:m]
loglik[1] = sum(logmarglik(vec(mcmcobj.n1s[cl[1], :]), mcmcobj.v[cl[1]],
                           priorC.A[linearidx], priorC.B[linearidx]))

linearidx = [(mcmcobj.C[cl[2], b] + 4 * (b - 1))::Int for b in 1:m]
loglik[2] = sum(logmarglik(vec(mcmcobj.n1s[cl[2], :]), mcmcobj.v[cl[2]],
                           priorC.A[linearidx], priorC.B[linearidx]))

ll = loglik[1] + loglik[2]
@test_approx_eq_eps mcmcobj.loglik ll ε
