# This file is part of Kpax3. License is MIT.

ifile = "data/proper_aa.fasta"
ofile = "../build/test.bin"

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

(m, n) = size(data)

# test constructor
settings = KSettings(ifile, ofile, maxclust=1, maxunit=1)

R = [13; 13; 42; 42; 76; 76]

k = length(unique(R))
maxclust = max(k, min(n, settings.maxclust))

emptycluster = trues(maxclust)
emptycluster[1:k] = false

cl = zeros(Int, maxclust)
cl[1:k] = find(!emptycluster)

v = zeros(Int, maxclust)
v[1:k] = [2; 2; 2]

n1s = zeros(Float64, maxclust, m)
n1s[1, :] = vec(sum(float(data[:, R .== 13]), 2))
n1s[2, :] = vec(sum(float(data[:, R .== 42]), 2))
n1s[3, :] = vec(sum(float(data[:, R .== 76]), 2))

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, settings.γ, settings.r,
                           maxclust=settings.maxclust)

# priorC should be resized, test that as well
@test length(priorC.logω) == 1

state = AminoAcidState(data, R, priorR, priorC, settings)

@test length(priorC.logω) == k

@test state.R == [1; 1; 2; 2; 3; 3]

@test isa(state.C, Matrix{UInt8})
@test size(state.C, 1) == maxclust
@test size(state.C, 2) == m

@test state.emptycluster == emptycluster
@test state.cl == cl
@test state.k == k

@test state.v == v
@test state.n1s == n1s
@test state.unit == Vector{Int}[sum(state.R .== g) > 0 ?
                                find(state.R .== g) :
                                [0] for g in 1:maxclust]

@test state.logpR == logdPriorRow(n, k, v, priorR)
@test_approx_eq_eps state.logpC[1] logpriorC(state.C, state.cl, state.k,
                                             priorC) ε
@test_approx_eq_eps state.logpC[2] logcondpostC(state.C, state.cl,
                                                state.k, state.v,
                                                state.n1s, priorC) ε

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

@test_approx_eq_eps state.logpp (logdPriorRow(n, k, v, priorR) +
                                 logpriorC(state.C, state.cl, state.k, priorC) +
                                 ll) ε

# test methods
resizelogω!(priorC, 4)

# test resizestate!
settings = KSettings(ifile, ofile, maxclust=1, maxunit=1)
state1 = AminoAcidState(data, [1; 1; 1; 1; 1; 1], priorR, priorC, settings)

len = 3

R = copy(state1.R)
k = state1.k

C = zeros(UInt8, len, m)
C[1, :] = copy(vec(state1.C))

emptycluster = trues(len)
emptycluster[1] = false

cl = zeros(Int, len)
cl[1] = 1

v = zeros(Int, len)
v[1] = 6

n1s = zeros(Float64, len, m)
n1s[1, :] = copy(vec(state1.n1s))

unit = Vector{Int}[sum(state1.R .== g) > 0 ? find(state1.R .== g) :
                   zeros(Int, n) for g in 1:len]

logpR = state1.logpR
logpC = copy(state1.logpC)
loglik = state1.loglik
logpp = state1.logpp

resizestate!(state1, 3)

@test state1.R == R
@test state1.C == C
@test state1.emptycluster == emptycluster
@test state1.cl == cl
@test state1.k == k
@test state1.v == v
@test state1.n1s == n1s
@test state1.unit == unit
@test state1.logpR == logpR
@test state1.logpC == logpC
@test state1.loglik == loglik
@test state1.logpp == logpp

settings = KSettings(ifile, ofile, maxclust=2, maxunit=1)
state1 = AminoAcidState(data, [1; 1; 1; 1; 1; 1], priorR, priorC, settings)

len = 4

R = copy(state1.R)
k = state1.k

C = zeros(UInt8, len, m)
C[1, :] = copy(vec(state1.C[1, :]))

emptycluster = trues(len)
emptycluster[1] = false

cl = zeros(Int, len)
cl[1] = 1

v = zeros(Int, len)
v[1] = 6

n1s = zeros(Float64, len, m)
n1s[1, :] = copy(vec(state1.n1s[1, :]))

unit = Vector{Int}[sum(state1.R .== g) > 0 ? find(state1.R .== g) : [0]
                   for g in 1:len]

logpR = state1.logpR
logpC = copy(state1.logpC)
loglik = state1.loglik
logpp = state1.logpp

resizestate!(state1, 3, settings)

@test state1.R == R
@test state1.C == C
@test state1.emptycluster == emptycluster
@test state1.cl == cl
@test state1.k == k
@test state1.v == v
@test state1.n1s == n1s
@test state1.unit == unit
@test state1.logpR == logpR
@test state1.logpC == logpC
@test state1.loglik == loglik
@test state1.logpp == logpp

# test copy functions
R = [3; 3; 1; 1; 5; 5]
g = sort(unique(R))
k = length(g)
len = maximum(g)

C = ones(UInt8, len, m)

emptycluster = trues(len)
emptycluster[g] = false

cl = zeros(Int, len)
cl[1:k] = g

v = zeros(Int, len)
v[g] = [2; 2; 2]

n1s = zeros(Float64, len, m)
n1s[1, :] = vec(sum(float(data[:, R .== 1]), 2))
n1s[3, :] = vec(sum(float(data[:, R .== 3]), 2))
n1s[5, :] = vec(sum(float(data[:, R .== 5]), 2))

unit = Vector{Int}[sum(R .== g) > 0 ? find(R .== g) : [0] for g in 1:n]

logpR = -6.5792512120101012129680384532548487186431884765625
logpC = [-9.1948612277878307708078864379785954952239990234375;
         -7.2468962917275181467857692041434347629547119140625]
loglik = -66.7559125873377894322402426041662693023681640625
logpp = logpR + logpC[1] + loglik

# test basic copying
state1 = AminoAcidState(copy(R), copy(C), copy(emptycluster), copy(cl), copy(k),
                        copy(v), copy(n1s), deepcopy(unit), logpR, copy(logpC),
                        loglik, logpp)

state2 = copystate(state1)

@test state2.R == R
@test state2.C == C
@test state2.emptycluster == emptycluster
@test state2.cl == cl
@test state2.k == k
@test state2.v == v
@test state2.n1s == n1s
@test state2.unit == unit
@test state2.logpR == logpR
@test state2.logpC == logpC
@test state2.loglik == loglik
@test state2.logpp == logpp

# objects should be copies, not references to state1 objects
# change state1 objects and see if state2 changes as well
state1.R = [1; 1; 1; 1; 2; 3]
fill!(state1.C, UInt8(2))
state1.emptycluster = [false; false; false; true; true]
state1.cl = [1; 2; 3; 0; 0]
state1.k = 3
state1.v = [4; 1; 1; 0; 2]
state1.n1s[1, :] = vec(sum(float(data[:, state1.R .== 1]), 2))
state1.n1s[2, :] = vec(sum(float(data[:, state1.R .== 2]), 2))
state1.n1s[3, :] = vec(sum(float(data[:, state1.R .== 3]), 2))
state1.unit[1] = [1; 2; 3; 4]
state1.unit[2] = [5; 0]
state1.unit[3] = [6]
state1.logpR = logdPriorRow(n, state1.k, state1.v[1:state1.k], priorR)
state1.logpC[1] = logpriorC(state1.C, state1.cl, state1.k, priorC)
state1.logpC[2] = logcondpostC(state1.C, state1.cl, state1.k, state1.v,
                               state1.n1s, priorC)
state1.loglik = loglikelihood(state1.C, state1.cl, state1.k, state1.v,
                              state1.n1s, priorC)
state1.logpp = state1.logpR + state1.logpC[1] + state1.loglik

@test state2.R == R
@test state2.C == C
@test state2.emptycluster == emptycluster
@test state2.cl == cl
@test state2.k == k
@test state2.v == v
@test state2.n1s == n1s
@test state2.unit == unit
@test state2.logpR == logpR
@test state2.logpC == logpC
@test state2.loglik == loglik
@test state2.logpp == logpp

# test copying in place
# state2 should be resized
R = [3; 3; 1; 1; 5; 5]
g = sort(unique(R))
k = length(g)
len = maximum(g)

C = ones(UInt8, len, m)

emptycluster = trues(len)
emptycluster[g] = false

cl = zeros(Int, len)
cl[1:k] = g

v = zeros(Int, len)
v[g] = [2; 2; 2]

n1s = zeros(Float64, len, m)
n1s[1, :] = vec(sum(float(data[:, R .== 1]), 2))
n1s[3, :] = vec(sum(float(data[:, R .== 3]), 2))
n1s[5, :] = vec(sum(float(data[:, R .== 5]), 2))

unit = Vector{Int}[sum(R .== g) > 0 ? find(R .== g) : [0] for g in 1:n]

logpR = -6.5792512120101012129680384532548487186431884765625
logpC = [-9.1948612277878307708078864379785954952239990234375;
         -7.2468962917275181467857692041434347629547119140625]
loglik = -66.7559125873377894322402426041662693023681640625
logpp = logpR + logpC[1] + loglik

# test basic copying
state1 = AminoAcidState(copy(R), copy(C), copy(emptycluster), copy(cl), copy(k),
                        copy(v), copy(n1s), deepcopy(unit), logpR, copy(logpC),
                        loglik, logpp)

settings = KSettings(ifile, ofile, maxclust=1, maxunit=1)
state2 = AminoAcidState(data, [1; 1; 1; 1; 1; 1], priorR, priorC, settings)

copystate!(state2, state1)

l = state1.cl[1:state1.k]

@test state2.R == state1.R
@test state2.C[l, :] == state1.C[l, :]
@test state2.emptycluster[l] == state1.emptycluster[l]
@test state2.cl[1:state2.k] == state1.cl[1:state1.k]
@test state2.k == state1.k
@test state2.v[l] == state1.v[l]
@test state2.n1s[l, :] == state1.n1s[l, :]
for g in l
  @test state2.unit[g][1:state2.v[g]] == state1.unit[g][1:state1.v[g]]
end
@test state2.logpR == state1.logpR
@test state2.logpC == state1.logpC
@test state2.loglik == state1.loglik
@test state2.logpp == state1.logpp

state1.R = [1; 3; 3; 4; 4; 5]
fill!(state1.C, UInt8(2))
state1.emptycluster = [false; true; false; false; true; false]
state1.cl = [1; 3; 4; 5; 0; 0]
state1.k = 4
state1.v = [1; 0; 2; 2; 1; 0]
fill!(state1.n1s, 0.0)
state1.n1s[1, :] = vec(sum(float(data[:, state1.R .== 1]), 2))
state1.n1s[3, :] = vec(sum(float(data[:, state1.R .== 3]), 2))
state1.n1s[4, :] = vec(sum(float(data[:, state1.R .== 4]), 2))
state1.n1s[5, :] = vec(sum(float(data[:, state1.R .== 5]), 2))
state1.unit[1] = [1; 0; 0; 0; 0; 0]
state1.unit[2] = [0; 0; 0; 0; 0; 0]
state1.unit[3] = [2; 3; 0; 0; 0; 0]
state1.unit[4] = [4; 5; 0; 0; 0; 0]
state1.unit[5] = [6; 0; 0; 0; 0; 0]
state1.logpR = logdPriorRow(n, state1.k, state1.v[[1; 3; 4; 5]], priorR)
state1.logpC[1] = logpriorC(state1.C, state1.cl, state1.k, priorC)
state1.logpC[2] = logcondpostC(state1.C, state1.cl, state1.k, state1.v,
                               state1.n1s, priorC)
state1.loglik = loglikelihood(state1.C, state1.cl, state1.k, state1.v,
                              state1.n1s, priorC)
state1.logpp = state1.logpR + state1.logpC[1] + state1.loglik

@test state2.R != state1.R
@test state2.C[l, :] != state1.C[l, :]
@test state2.emptycluster[l] != state1.emptycluster[l]
@test state2.cl[1:state2.k] != state1.cl[1:state2.k]
@test state2.k != state1.k
@test state2.v[l] != state1.v[l]
@test state2.n1s[l, :] != state1.n1s[l, :]
for g in l
  @test state2.unit[g][1:state2.v[g]] != state1.unit[g][1:state1.v[g]]
end
@test state2.logpR != state1.logpR
@test state2.logpC != state1.logpC
@test state2.loglik != state1.loglik
@test state2.logpp != state1.logpp

# state2 should be fine without resize
R = [3; 3; 1; 1; 5; 5]
g = sort(unique(R))
k = length(g)
len = maximum(g)

C = ones(UInt8, len, m)

emptycluster = trues(len)
emptycluster[g] = false

cl = zeros(Int, len)
cl[1:k] = g

v = zeros(Int, len)
v[g] = [2; 2; 2]

n1s = zeros(Float64, len, m)
n1s[1, :] = vec(sum(float(data[:, R .== 1]), 2))
n1s[3, :] = vec(sum(float(data[:, R .== 3]), 2))
n1s[5, :] = vec(sum(float(data[:, R .== 5]), 2))

unit = Vector{Int}[sum(R .== g) > 0 ? find(R .== g) : [0] for g in 1:n]

logpR = -6.5792512120101012129680384532548487186431884765625
logpC = [-9.1948612277878307708078864379785954952239990234375;
         -7.2468962917275181467857692041434347629547119140625]
loglik = -66.7559125873377894322402426041662693023681640625
logpp = logpR + logpC[1] + loglik

# test basic copying
state1 = AminoAcidState(copy(R), copy(C), copy(emptycluster), copy(cl), copy(k),
                        copy(v), copy(n1s), deepcopy(unit), logpR, copy(logpC),
                        loglik, logpp)

settings = KSettings(ifile, ofile, maxclust=6, maxunit=6)
state2 = AminoAcidState(data, [1; 1; 1; 1; 1; 1], priorR, priorC, settings)

copystate!(state2, state1)

l = state1.cl[1:state1.k]

@test state2.R == state1.R
@test state2.C[l, :] == state1.C[l, :]
@test state2.emptycluster[l] == state1.emptycluster[l]
@test state2.cl[1:state2.k] == state1.cl[1:state1.k]
@test state2.k == state1.k
@test state2.v[l] == state1.v[l]
@test state2.n1s[l, :] == state1.n1s[l, :]
for g in l
  @test state2.unit[g][1:state2.v[g]] == state1.unit[g][1:state1.v[g]]
end
@test state2.logpR == state1.logpR
@test state2.logpC == state1.logpC
@test state2.loglik == state1.loglik
@test state2.logpp == state1.logpp

state1.R = [1; 3; 3; 4; 4; 5]
fill!(state1.C, UInt8(2))
state1.emptycluster = [false; true; false; false; true; false]
state1.cl = [1; 3; 4; 5; 0; 0]
state1.k = 4
state1.v = [1; 0; 2; 2; 1; 0]
fill!(state1.n1s, 0.0)
state1.n1s[1, :] = vec(sum(float(data[:, state1.R .== 1]), 2))
state1.n1s[3, :] = vec(sum(float(data[:, state1.R .== 3]), 2))
state1.n1s[4, :] = vec(sum(float(data[:, state1.R .== 4]), 2))
state1.n1s[5, :] = vec(sum(float(data[:, state1.R .== 5]), 2))
state1.unit[1] = [1; 0; 0; 0; 0; 0]
state1.unit[2] = [0; 0; 0; 0; 0; 0]
state1.unit[3] = [2; 3; 0; 0; 0; 0]
state1.unit[4] = [4; 5; 0; 0; 0; 0]
state1.unit[5] = [6; 0; 0; 0; 0; 0]
state1.logpR = logdPriorRow(n, state1.k, state1.v[[1; 3; 4; 5]], priorR)
state1.logpC[1] = logpriorC(state1.C, state1.cl, state1.k, priorC)
state1.logpC[2] = logcondpostC(state1.C, state1.cl, state1.k, state1.v,
                               state1.n1s, priorC)
state1.loglik = loglikelihood(state1.C, state1.cl, state1.k, state1.v,
                              state1.n1s, priorC)
state1.logpp = state1.logpR + state1.logpC[1] + state1.loglik

@test state2.R != state1.R
@test state2.C[l, :] != state1.C[l, :]
@test state2.emptycluster[l] != state1.emptycluster[l]
@test state2.cl[1:state2.k] != state1.cl[1:state2.k]
@test state2.k != state1.k
@test state2.v[l] != state1.v[l]
@test state2.n1s[l, :] != state1.n1s[l, :]
for g in l
  @test state2.unit[g][1:state2.v[g]] != state1.unit[g][1:state1.v[g]]
end
@test state2.logpR != state1.logpR
@test state2.logpC != state1.logpC
@test state2.loglik != state1.loglik
@test state2.logpp != state1.logpp

# test updating in place
# state2 should be resized
settings = KSettings(ifile, ofile, maxclust=1, maxunit=1)
state1 = AminoAcidState(data, [1; 1; 1; 1; 2; 3], priorR, priorC, settings)
state2 = AminoAcidState(data, [1; 1; 1; 1; 1; 1], priorR, priorC, settings)

updatestate!(state2, data, [1; 1; 1; 1; 2; 3], priorR, priorC, settings)

@test state2.R == state1.R
@test state2.C == state1.C
@test state2.emptycluster == state1.emptycluster
@test state2.cl == state1.cl
@test state2.k == state1.k
@test state2.v == state1.v
@test state2.n1s == state1.n1s
@test state2.unit[1] == [1; 2; 3; 4; 5; 6]
@test state2.unit[2] == [5]
@test state2.unit[3] == [6]
@test state2.logpR == state1.logpR
@test state2.logpC == state1.logpC
@test state2.loglik == state1.loglik
@test state2.logpp == state1.logpp

settings = KSettings(ifile, ofile, maxclust=6, maxunit=6)
state1 = AminoAcidState(data, [1; 1; 1; 2; 2; 2], priorR, priorC, settings)
state2 = AminoAcidState(data, [1; 1; 2; 2; 3; 3], priorR, priorC, settings)

updatestate!(state2, data, [1; 1; 1; 2; 2; 2], priorR, priorC, settings)

l = state1.cl[1:state1.k]

@test state2.R == state1.R
@test state2.C[l, :] == state1.C[l, :]
@test state2.emptycluster == state1.emptycluster
@test state2.cl[l] == state1.cl[l]
@test state2.k == state1.k
@test state2.v[l] == state1.v[l]
@test state2.n1s[l, :] == state1.n1s[l, :]
@test state2.unit[1] == [1; 2; 3; 0; 0; 0]
@test state2.unit[2] == [4; 5; 6; 0; 0; 0]
@test state2.unit[3] == [5; 6; 0; 0; 0; 0]
@test state2.logpR == state1.logpR
@test state2.logpC == state1.logpC
@test state2.loglik == state1.loglik
@test state2.logpp == state1.logpp

# TODO: How to test properly initializestate?
settings = KSettings(ifile, ofile, maxclust=6, maxunit=6)

D = zeros(Float64, n, n)
for j in 1:(n - 1), i in (j + 1):n
  D[i, j] = D[j, i] = sum(data[:, j] .!= data[:, i]) / m
end

s = initializestate(data, D, 1:6, priorR, priorC, settings)

@test isa(s.R, Vector{Int})
@test all(s.R .> 0)

t = AminoAcidState(data, s.R, priorR, priorC, settings)

l = t.cl[1:t.k]

@test s.R == t.R
@test s.C == t.C
@test s.emptycluster == t.emptycluster
@test s.cl == t.cl
@test s.k == t.k
@test s.v == t.v
@test s.n1s == t.n1s
for g in l
  @test s.unit[g][1:s.v[g]] == t.unit[g][1:t.v[g]]
end
@test s.logpR == t.logpR
@test s.logpC == t.logpC
@test s.loglik == t.loglik

@test_approx_eq_eps s.logpp t.logpp ε
