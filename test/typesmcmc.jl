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

R = [1; 1; 1; 1; 2; 2]
v = [4; 2]
n1s = [0.0 1.0;
       4.0 1.0;
       1.0 2.0;
       3.0 0.0;
       2.0 0.0;
       1.0 1.0;
       3.0 0.0;
       1.0 2.0;
       1.0 0.0;
       2.0 1.0;
       1.0 1.0;
       0.0 1.0;
       3.0 0.0;
       1.0 1.0;
       2.0 2.0;
       2.0 0.0;
       3.0 0.0;
       1.0 2.0]

k = maximum(R)

α = 0.0
θ = 1.0
priorR = EwensPitman(α, θ)

g = [0.6; 0.35; 0.05]
r = log(0.001) / log(0.95)
priorC = AminoAcidPriorCol(data, k, g, r)

settings = KSettings("typesmcmc.bin", 1, 0, 1, [1.0; 0.0; 0.0], 0.0, 1.0,
                     [0.6; 0.35; 0.05], 135.0, 1.0, 1.0, 5.0, 3, 1, true, 1)

mcmcobj = AminoAcidMCMC(data, R, priorR, priorC, settings)

@test mcmcobj.R == R

@test isa(mcmcobj.C, Array{UInt8, 2})
@test size(mcmcobj.C, 1) == settings.maxclust
@test size(mcmcobj.C, 2) == m

@test mcmcobj.filledcluster == [true; true; false]
@test mcmcobj.cl == [1; 2]

@test mcmcobj.v == [v; 0]
@test mcmcobj.n1s == vcat(n1s', zeros(Float64, 1, m))
@test mcmcobj.unit == Vector{Int}[sum(R .== g) > 0 ? find(R .== g) : [0]
                                  for g in 1:3]

@test mcmcobj.logprR == logdPriorRow(n, k, v, priorR)
@test_approx_eq_eps mcmcobj.logprC logpriorC(mcmcobj.C, mcmcobj.cl, priorC.logγ,
                                             priorC.logω) eps()
@test_approx_eq_eps mcmcobj.logpocC logcondpostC(mcmcobj.C, mcmcobj.cl,
                                                 mcmcobj.v, mcmcobj.n1s,
                                                 priorC.logγ, priorC.logω,
                                                 priorC.A, priorC.B) eps()

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
loglik[1] = sum(logmarglik(n1s[:, 1], v[1], priorC.A[linearidx],
                           priorC.B[linearidx]))

linearidx = [(C[2, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[2] = sum(logmarglik(n1s[:, 2], v[2], priorC.A[linearidx],
                           priorC.B[linearidx]))

@test loglik[1] == sum(logmarglik(n1s[:, 1], v[1], A[:, 1], B[:, 1]))
@test loglik[2] == sum(logmarglik(n1s[:, 2], v[2], A[:, 2], B[:, 2]))

loglik = zeros(Float64, 2)
linearidx = [(mcmcobj.C[1, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[1] = sum(logmarglik(vec(mcmcobj.n1s[1, :]), mcmcobj.v[1],
                           priorC.A[linearidx], priorC.B[linearidx]))

linearidx = [(mcmcobj.C[2, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[2] = sum(logmarglik(vec(mcmcobj.n1s[2, :]), mcmcobj.v[2],
                           priorC.A[linearidx], priorC.B[linearidx]))

ll = loglik[1] + loglik[2]
@test_approx_eq_eps mcmcobj.loglik ll ε
