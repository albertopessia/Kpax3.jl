# This file is part of Kpax3. License is MIT.

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
                     [0.6; 0.35; 0.05], 135.0, 1.0, 1.0, 5.0, 2, 1, true, 1)

mcmcobj = AminoAcidMCMC(data, R, priorR, priorC, settings)

@test mcmcobj.R == R

@test isa(mcmcobj.C, Array{UInt8, 2})
@test size(mcmcobj.C, 1) == settings.maxclust
@test size(mcmcobj.C, 2) == m

@test isa(mcmcobj.cluster, Array{KCluster, 1})
@test length(mcmcobj.cluster) == settings.maxclust

linearidx = [(mcmcobj.C[1, b] + 4 * (b - 1))::Int for b in 1:m]

@test mcmcobj.cluster[1].v == v[1]
@test mcmcobj.cluster[1].unit == [1; 2; 3; 4; zeros(Int, settings.maxclust - 4)]
@test mcmcobj.cluster[1].n1s == n1s[:, 1]
@test mcmcobj.cluster[1].loglik == sum(logmarglik(n1s[:, 1], v[1],
                                       priorC.A[linearidx],
                                       priorC.B[linearidx]))

linearidx = [(mcmcobj.C[2, b] + 4 * (b - 1))::Int for b in 1:m]

@test mcmcobj.cluster[2].v == v[2]
@test mcmcobj.cluster[2].unit == [5; 6; zeros(Int, settings.maxclust - 2)]
@test mcmcobj.cluster[2].n1s == n1s[:, 2]
@test mcmcobj.cluster[2].loglik == sum(logmarglik(n1s[:, 2], v[2],
                                       priorC.A[linearidx],
                                       priorC.B[linearidx]))

@test mcmcobj.emptycluster == [false; false; trues(settings.maxclust - 2)]

@test mcmcobj.k == k

@test mcmcobj.logprR == logdPriorRow(n, k, v, priorR)
@test mcmcobj.logprC == logpriorC(mcmcobj.C, mcmcobj.emptycluster,
                                  priorC.logγ, priorC.logω)

@test mcmcobj.loglik == mcmcobj.cluster[1].loglik + mcmcobj.cluster[2].loglik

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
