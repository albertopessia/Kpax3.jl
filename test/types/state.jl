# This file is part of Kpax3. License is MIT.

data = UInt8[0 0 0 0 0 1;
             1 1 1 1 1 0;
             0 0 1 0 1 1;
             1 1 0 1 0 0;
             1 1 0 0 0 0;
             0 0 0 1 1 0;
             1 1 1 0 0 0;
             0 0 0 1 1 1;
             0 0 1 0 0 0;
             1 0 0 1 0 1;
             0 1 0 0 1 0;
             0 0 0 0 0 1;
             1 1 1 0 0 0;
             0 0 0 1 1 0;
             1 1 0 0 1 1;
             0 0 1 1 0 0;
             1 1 0 1 0 0;
             0 0 1 0 1 1]

m, n = size(data)

settings = KSettings("../build/test.bin", T=1, burnin=0, tstep=1,
                     op=[0.6; 0.3; 0.1], α=0.0, θ=1.0, γ=[0.6; 0.35; 0.05],
                     r=135.0, λs1=1.0, λs2=1.0, parawm=5.0, maxclust=100,
                     maxunit=1, verbose=true, verbosestep=1)

R = [13; 13; 42; 42; 76; 76]

k = length(unique(R))

emptycluster = trues(n)
emptycluster[1:k] = false

cl = zeros(Int, n)
cl[1:k] = find(!emptycluster)

v = zeros(Int, n)
v[1:k] = [2; 2; 2]

n1s = zeros(Float64, n, m)
n1s[1, :] = vec(sum(float(data[:, R.== 13]), 2))
n1s[2, :] = vec(sum(float(data[:, R.== 42]), 2))
n1s[3, :] = vec(sum(float(data[:, R.== 76]), 2))

α = 0.0
θ = 1.0
priorR = EwensPitman(α, θ)

g = [1.0; 1.0; 1.0]
r = log(0.001) / log(0.95)
priorC = AminoAcidPriorCol(data, k, g, r)

state = AminoAcidState(data, R, priorR, priorC, settings)

@test state.R == [1; 1; 2; 2; 3; 3]

@test isa(state.C, Array{UInt8, 2})
@test size(state.C, 1) == n
@test size(state.C, 2) == m

@test state.emptycluster == emptycluster
@test state.cl == cl
@test state.k == k

@test state.v == v
@test state.n1s == n1s
@test state.unit == Vector{Int}[sum(state.R .== g) > 0 ?
                                find(state.R .== g) :
                                [0] for g in 1:n]

@test state.logpR == logdPriorRow(n, k, v, priorR)
@test_approx_eq_eps state.logpC[1] logpriorC(state.C, state.cl, state.k,
                                             priorC.logγ, priorC.logω) ε
@test_approx_eq_eps state.logpC[2] logcondpostC(state.C, state.cl,
                                                state.k, state.v,
                                                state.n1s, priorC.logω,
                                                priorC) ε

# is linearidx approach correct?
C = UInt8[1 1 1 4 2 1 3 3 1 1 1 1 1 1 1 1 1 2;
          1 1 1 3 2 1 3 4 1 1 1 1 1 1 1 1 1 2;
          1 1 1 3 2 1 4 3 1 1 1 1 1 1 1 1 1 2]

A = [priorC.A[1,  1] priorC.A[1,  1] priorC.A[1,  1];
     priorC.A[1,  2] priorC.A[1,  2] priorC.A[1,  2];
     priorC.A[1,  3] priorC.A[1,  3] priorC.A[1,  3];
     priorC.A[4,  4] priorC.A[3,  4] priorC.A[3,  4];
     priorC.A[2,  5] priorC.A[2,  5] priorC.A[2,  5];
     priorC.A[1,  6] priorC.A[1,  6] priorC.A[1,  6];
     priorC.A[3,  7] priorC.A[3,  7] priorC.A[4,  7];
     priorC.A[3,  8] priorC.A[4,  8] priorC.A[3,  8];
     priorC.A[1,  9] priorC.A[1,  9] priorC.A[1,  9];
     priorC.A[1, 10] priorC.A[1, 10] priorC.A[1, 10];
     priorC.A[1, 11] priorC.A[1, 11] priorC.A[1, 11];
     priorC.A[1, 12] priorC.A[1, 12] priorC.A[1, 12];
     priorC.A[1, 13] priorC.A[1, 13] priorC.A[1, 13];
     priorC.A[1, 14] priorC.A[1, 14] priorC.A[1, 14];
     priorC.A[1, 15] priorC.A[1, 15] priorC.A[1, 15];
     priorC.A[1, 16] priorC.A[1, 16] priorC.A[1, 16];
     priorC.A[1, 17] priorC.A[1, 17] priorC.A[1, 17];
     priorC.A[2, 18] priorC.A[2, 18] priorC.A[2, 18]]

B = [priorC.B[1,  1] priorC.B[1,  1] priorC.B[1,  1];
     priorC.B[1,  2] priorC.B[1,  2] priorC.B[1,  2];
     priorC.B[1,  3] priorC.B[1,  3] priorC.B[1,  3];
     priorC.B[4,  4] priorC.B[3,  4] priorC.B[3,  4];
     priorC.B[2,  5] priorC.B[2,  5] priorC.B[2,  5];
     priorC.B[1,  6] priorC.B[1,  6] priorC.B[1,  6];
     priorC.B[3,  7] priorC.B[3,  7] priorC.B[4,  7];
     priorC.B[3,  8] priorC.B[4,  8] priorC.B[3,  8];
     priorC.B[1,  9] priorC.B[1,  9] priorC.B[1,  9];
     priorC.B[1, 10] priorC.B[1, 10] priorC.B[1, 10];
     priorC.B[1, 11] priorC.B[1, 11] priorC.B[1, 11];
     priorC.B[1, 12] priorC.B[1, 12] priorC.B[1, 12];
     priorC.B[1, 13] priorC.B[1, 13] priorC.B[1, 13];
     priorC.B[1, 14] priorC.B[1, 14] priorC.B[1, 14];
     priorC.B[1, 15] priorC.B[1, 15] priorC.B[1, 15];
     priorC.B[1, 16] priorC.B[1, 16] priorC.B[1, 16];
     priorC.B[1, 17] priorC.B[1, 17] priorC.B[1, 17];
     priorC.B[2, 18] priorC.B[2, 18] priorC.B[2, 18]]

loglik = zeros(Float64, 3)

linearidx = [(C[1, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[1] = sum(logmarglik(vec(n1s[cl[1], :]), v[cl[1]], priorC.A[linearidx],
                           priorC.B[linearidx]))

linearidx = [(C[2, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[2] = sum(logmarglik(vec(n1s[cl[2], :]), v[cl[2]], priorC.A[linearidx],
                           priorC.B[linearidx]))

linearidx = [(C[3, b] + 4 * (b - 1))::Int for b in 1:m]
loglik[3] = sum(logmarglik(vec(n1s[cl[3], :]), v[cl[3]], priorC.A[linearidx],
                           priorC.B[linearidx]))

@test loglik[1] == sum(logmarglik(vec(n1s[cl[1], :]), v[cl[1]], A[:, 1],
                                  B[:, 1]))
@test loglik[2] == sum(logmarglik(vec(n1s[cl[2], :]), v[cl[2]], A[:, 2],
                                  B[:, 2]))
@test loglik[3] == sum(logmarglik(vec(n1s[cl[3], :]), v[cl[3]], A[:, 3],
                                  B[:, 3]))

loglik = zeros(Float64, 3)

linearidx = [(state.C[cl[1], b] + 4 * (b - 1))::Int for b in 1:m]
loglik[1] = sum(logmarglik(vec(state.n1s[cl[1], :]), state.v[cl[1]],
                           priorC.A[linearidx], priorC.B[linearidx]))

linearidx = [(state.C[cl[2], b] + 4 * (b - 1))::Int for b in 1:m]
loglik[2] = sum(logmarglik(vec(state.n1s[cl[2], :]), state.v[cl[2]],
                           priorC.A[linearidx], priorC.B[linearidx]))

linearidx = [(state.C[cl[3], b] + 4 * (b - 1))::Int for b in 1:m]
loglik[3] = sum(logmarglik(vec(state.n1s[cl[3], :]), state.v[cl[3]],
                           priorC.A[linearidx], priorC.B[linearidx]))

ll = loglik[1] + loglik[2] + loglik[3]
@test_approx_eq_eps state.loglik ll ε
