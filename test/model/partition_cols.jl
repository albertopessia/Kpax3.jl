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

(m, n) = size(data)
n1s = Float64[1; 5; 3; 3; 2; 2; 3; 3; 1; 3; 2; 1; 3; 2; 4; 2; 3; 3]

r1 = 2.0
r2 = 100.0

A1 = zeros(Float64, 4, m)
A1[1, :] = (r1 + 1.0) * (n1s + 0.5) / (n + 1)
A1[2, :] = 1.0
A1[3, :] = 1.0
A1[4, :] = r1

B1 = zeros(Float64, 4, m)
B1[1, :] = (r1 + 1.0) - A1[1, :]
B1[2, :] = 1.0
B1[3, :] = r1
B1[4, :] = 1.0

A2 = zeros(Float64, 4, m)
A2[1, :] = n1s + 0.5
A2[2, :] = 1.0
A2[3, :] = 1.0
A2[4, :] = r2

B2 = zeros(Float64, 4, m)
B2[1, :] = n - n1s + 0.5
B2[2, :] = 1.0
B2[3, :] = r2
B2[4, :] = 1.0

for k in 1:n
  for γ in ([1.0; 0.0; 0.0], [0.0; 1.0; 0.0], [0.0; 0.0; 1.0],
            [0.4; 0.3; 0.3], [0.5; 0.3; 0.2], [0.7; 0.2; 0.1],
            [0.1; 0.1; 0.1], [0.3; 0.1; 0.1], [0.0; 0.2; 0.1])
    x1 = AminoAcidPriorCol(data, γ, r1, maxclust=k)
    x2 = AminoAcidPriorCol(data, γ, r2, maxclust=k)

    @test_approx_eq_eps x1.logγ[1] log(γ[1] / sum(γ)) ε
    @test_approx_eq_eps x2.logγ[2] log(γ[2] / sum(γ)) ε
    @test_approx_eq_eps x2.logγ[3] log(γ[3] / sum(γ)) ε

    @test_approx_eq_eps x1.logω[k][1] log(1.0 - 1.0 / k) ε
    @test_approx_eq_eps x2.logω[k][2] log(1.0 / k) ε

    for b in 1:m, s in 1:4
      @test_approx_eq_eps x1.A[s, b] A1[s, b] ε
      @test_approx_eq_eps x1.B[s, b] B1[s, b] ε

      @test_approx_eq_eps x2.A[s, b] A2[s, b] ε
      @test_approx_eq_eps x2.B[s, b] B2[s, b] ε
    end
  end
end

ifile = "data/proper_aa.fasta"
ofile = "../build/test.bin"

settings = KSettings(ifile, ofile, maxclust=100, maxunit=1)

data = UInt8[1 1 0 0 1 1;
             0 0 1 1 0 0;
             1 1 0 0 0 1;
             1 0 0 0 1 1]

(m, n) = size(data)

R = [1; 1; 2; 2; 3; 1]
k = length(unique(R))

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, settings.γ, settings.r)

no = AminoAcidState(data, R, priorR, priorC, settings)

C = zeros(UInt8, settings.maxclust, m)

cs = (UInt8[1; 1; 1], UInt8[2; 2; 2], UInt8[3; 3; 3], UInt8[3; 3; 4],
      UInt8[3; 4; 3], UInt8[4; 3; 3], UInt8[3; 4; 4], UInt8[4; 3; 4],
      UInt8[4; 4; 3], UInt8[4; 4; 4])

logp1 = zeros(Float64, (2 + 2^k)^m)
logp2 = zeros(Float64, (2 + 2^k)^m)
logp3 = zeros(Float64, (2 + 2^k)^m)

l = 0

for c1 in cs, c2 in cs, c3 in cs, c4 in cs
  l += 1

  tmp = hcat(c1, c2, c3, c4)

  no.C[no.cl[1], :] = C[1, :] = tmp[1, :]
  no.C[no.cl[2], :] = C[2, :] = tmp[2, :]
  no.C[no.cl[3], :] = C[3, :] = tmp[3, :]

  logp1[l] = logpriorC(no.C, no.cl, k, priorC)
  logp2[l] = logpriorC(C, k, priorC)
  logp3[l] = logcondpostC(no.C, no.cl, k, no.v, no.n1s, priorC)
end

M = maximum(logp1)
p1 = exp(M + log(sum(exp(logp1 - M))))

M = maximum(logp2)
p2 = exp(M + log(sum(exp(logp2 - M))))

M = maximum(logp3)
p3 = exp(M + log(sum(exp(logp3 - M))))

@test_approx_eq_eps p1 1.0 ε
@test_approx_eq_eps p2 1.0 ε
@test_approx_eq_eps p3 1.0 ε

ss = UInt8[1; 2; 3]
logp4 = zeros(Float64, 3^m)
l = 0

for s1 in ss, s2 in ss, s3 in ss, s4 in ss
  l += 1
  S = [s1; s2; s3; s4]
  logp4[l] = logcondpostS(S, no.cl, k, no.v, no.n1s, priorC)
end

M = maximum(logp4)
p4 = exp(M + log(sum(exp(logp4 - M))))

@test_approx_eq_eps p4 1.0 ε

# suppose we merge cluster 1 and cluster 3
me = AminoAcidState(data, R, priorR, priorC, settings)

mergesupport = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))
mergevi = me.v[me.cl[3]] + me.v[me.cl[1]]
mergeni = vec(no.n1s[no.cl[3], :] + no.n1s[no.cl[1], :])

# suppose we split cluster 1
sp = AminoAcidState(data, R, priorR, priorC, settings)

splitsupport = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))
splitsupport.vi = sp.v[sp.cl[1]] - 1
splitsupport.ni = vec(sp.n1s[sp.cl[1], :]) - float(data[:, 2])
splitsupport.vj = 1
splitsupport.nj = float(data[:, 2])

# suppose we move unit 5 into cluster 1 (merge)
brw_1 = AminoAcidState(data, R, priorR, priorC, settings)
brwsupport_1 = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))

brwsupport_1.ni = float(data[:, 5])

# suppose we move unit 2 into cluster 3
brw_2 = AminoAcidState(data, R, priorR, priorC, settings)
brwsupport_2 = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))

brwsupport_2.ni = float(data[:, 6])

# suppose we move unit 2 into its own cluster (split)
brw_3 = AminoAcidState(data, R, priorR, priorC, settings)
brwsupport_3 = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))

brwsupport_3.ni = float(data[:, 2])

N = 1000000

Sp1 = zeros(Float64, 3, m)
Sp2 = zeros(Float64, 3, m)
Sp3 = zeros(Float64, 3, m)
Sp4 = zeros(Float64, 3, m)
Sp5 = zeros(Float64, 3, m)
Sp6 = zeros(Float64, 3, m)

Cp1 = zero(Float64)
Cp2 = zero(Float64)
Cp3 = zero(Float64)
Cp4 = zero(Float64)
Cp5 = zero(Float64)
Cp6 = zero(Float64)

Ctest1 = UInt8[2 1 4 1; 2 1 3 1; 2 1 3 1]
Ctest2 = UInt8[1 4 2 2; 1 3 2 2]
Ctest3 = UInt8[1 4 2 1; 1 3 2 1; 1 3 2 1; 1 3 2 1]
Ctest4 = UInt8[1 4 2 2; 1 3 2 2]
Ctest5 = UInt8[1 4 1 2; 1 3 1 2; 1 3 1 2]
Ctest6 = UInt8[1 4 2 1; 1 3 2 1; 1 3 2 1; 1 3 2 1]

#=
R = [1; 1; 2; 2; 3; 1]
k = length(unique(R))

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, settings.γ, settings.r)

ss = UInt8[1; 2; 3]
obj = AminoAcidState(data, R, priorR, priorC, settings)

M = -Inf
for s1 in ss, s2 in ss, s3 in ss, s4 in ss
  S = [s1; s2; s3; s4]
  lp = logcondpostS(S, obj.cl, obj.k, obj.v, obj.n1s, priorC)

  if lp > M
    M = lp
  end
end

trueSp = zeros(Float64, 3, 4)
for s1 in ss, s2 in ss, s3 in ss, s4 in ss
  S = [s1; s2; s3; s4]
  lp = logcondpostS(S, obj.cl, obj.k, obj.v, obj.n1s, priorC)

  trueSp[s1, 1] += exp(lp - M)
  trueSp[s2, 2] += exp(lp - M)
  trueSp[s3, 3] += exp(lp - M)
  trueSp[s4, 4] += exp(lp - M)
end

trueSp = exp(M + log(trueSp))

for i in 1:12
  @printf("%.100f\n", trueSp[i])
end

rpostpartitioncols!(obj.C, obj.cl, obj.k, obj.v, obj.n1s, priorC);
obj.C[obj.cl[1:obj.k], :]

obj.C[obj.cl[1:obj.k], :] = Ctest1

p = exp(logcondpostC(obj.C, obj.cl, obj.k, obj.v, obj.n1s, priorC))
@printf("%.100f\n", p)
=#

trueSp1 = hcat([0.514127006003589515081841909704962745308876037597656250000,
                0.389747596850238686716494385109399445354938507080078125000,
                0.096125397146171798201663705185637809336185455322265625000],
               [0.469929717191958307154209251166321337223052978515625000000,
                0.356242671218090556362056986472452990710735321044921875000,
                0.173827611589951080972582531103398650884628295898437500000],
               [0.400605930740164040138040491001447662711143493652343750000,
                0.402853999936236717438475807284703478217124938964843750000,
                0.196540069323599214667908086084935348480939865112304687500],
               [0.653732577496086042501133306359406560659408569335937500000,
                0.344352962631683023886353112175129354000091552734375000000,
                0.001914459872230926673619677558235707692801952362060546875])

trueSp2 = hcat([0.3886912881957179100034238672378705814480781555175781250000,
                0.4041020269672812581518428487470373511314392089843750000000,
                0.2072066848370019975789091404294595122337341308593750000000],
               [0.3886912881957179100034238672378705814480781555175781250000,
                0.4041020269672812581518428487470373511314392089843750000000,
                0.2072066848370019975789091404294595122337341308593750000000],
               [0.6316406983875068048561729483481030911207199096679687500000,
                0.3629630051340117513625216361106140539050102233886718750000,
                0.0053962964784825991071404160948077333159744739532470703125],
               [0.6316406983875068048561729483481030911207199096679687500000,
                0.3629630051340117513625216361106140539050102233886718750000,
                0.0053962964784825991071404160948077333159744739532470703125])

trueSp3 = hcat([0.624219233607703927191323600709438323974609375000000000000,
                0.354417607489966701717065689081209711730480194091796875000,
                0.021363158902329794364138848550283000804483890533447265625],
               [0.536479879971809014982397911808220669627189636230468750000,
                0.304601180625598710882684372336370870471000671386718750000,
                0.158918939402592662712976334660197608172893524169921875000],
               [0.508985055567922395347579822555417194962501525878906250000,
                0.417055884580573688058535708478302694857120513916015625000,
                0.073959059851504194149640625255415216088294982910156250000],
               [0.508985055567922284325277360039763152599334716796875000000,
                0.417055884580573688058535708478302694857120513916015625000,
                0.073959059851504277416367472142155747860670089721679687500])

trueSp4 = hcat([0.3886912881957179100034238672378705814480781555175781250000,
                0.4041020269672812581518428487470373511314392089843750000000,
                0.2072066848370019975789091404294595122337341308593750000000],
               [0.3886912881957179100034238672378705814480781555175781250000,
                0.4041020269672812581518428487470373511314392089843750000000,
                0.2072066848370019975789091404294595122337341308593750000000],
               [0.6316406983875068048561729483481030911207199096679687500000,
                0.3629630051340117513625216361106140539050102233886718750000,
                0.0053962964784825991071404160948077333159744739532470703125],
               [0.6316406983875068048561729483481030911207199096679687500000,
                0.3629630051340117513625216361106140539050102233886718750000,
                0.0053962964784825991071404160948077333159744739532470703125])

trueSp5 = hcat([0.525901931472022354796536092180758714675903320312500000000,
                0.372274609811358958566529508971143513917922973632812500000,
                0.101823458716617382124880464289162773638963699340820312500],
               [0.477315859308131362759297644515754655003547668457031250000,
                0.337881579524415998072583988687256351113319396972656250000,
                0.184802561167451473433942510382621549069881439208984375000],
               [0.612810148788300645961157897545490413904190063476562500000,
                0.382574359005774355946272180517553351819515228271484375000,
                0.004615492205923732611794196856180860777385532855987548828],
               [0.612810148788300645961157897545490413904190063476562500000,
                0.382574359005774355946272180517553351819515228271484375000,
                0.004615492205923739550688100763409238425083458423614501953])

trueSp6 = hcat([0.624219233607703927191323600709438323974609375000000000000,
                0.354417607489966701717065689081209711730480194091796875000,
                0.021363158902329794364138848550283000804483890533447265625],
               [0.536479879971809014982397911808220669627189636230468750000,
                0.304601180625598710882684372336370870471000671386718750000,
                0.158918939402592662712976334660197608172893524169921875000],
               [0.508985055567922395347579822555417194962501525878906250000,
                0.417055884580573688058535708478302694857120513916015625000,
                0.073959059851504194149640625255415216088294982910156250000],
               [0.508985055567922284325277360039763152599334716796875000000,
                0.417055884580573688058535708478302694857120513916015625000,
                0.073959059851504277416367472142155747860670089721679687500])

trueCp1 = 0.023444249569984754177909280770109035074710845947265625
trueCp2 = 0.0106092811116244335745140148219434195198118686676025390625
trueCp3 = 0.0209465448253591153549013625934094307012856006622314453125
trueCp4 = 0.0106092811116244335745140148219434195198118686676025390625
trueCp5 = 0.02277784562652716837671817984301014803349971771240234375
trueCp6 = 0.0209465448253591153549013625934094307012856006622314453125

for t in 1:N
  rpostpartitioncols!(no.C, no.cl, k, no.v, no.n1s, priorC)

  simcmerge!(k - 1, me.cl[3], me.cl[1], mergevi, mergeni, priorC,
             mergesupport, me)

  simcsplit!(k + 1, sp.cl[1], priorC, splitsupport, sp)

  simcbrw!(k - 1, 3, 1, priorC, brwsupport_1, brw_1)
  simcbrw!(k, 1, 3, priorC, brwsupport_2, brw_2)
  simcbrw!(k + 1, 1, 4, priorC, brwsupport_3, brw_3)

  for b in 1:m
    if no.C[no.cl[1], b] == 0x01
      @test no.C[no.cl[2], b] == 0x01
      @test no.C[no.cl[3], b] == 0x01

      Sp1[1, b] += 1.0
    elseif no.C[no.cl[1], b] == 0x02
      @test no.C[no.cl[2], b] == 0x02
      @test no.C[no.cl[3], b] == 0x02

      Sp1[2, b] += 1.0
    else
      Sp1[3, b] += 1.0
    end

    if mergesupport.C[1, b] == 0x01
      @test mergesupport.C[2, b] == 0x01

      Sp2[1, b] += 1.0
    elseif mergesupport.C[1, b] == 0x02
      @test mergesupport.C[2, b] == 0x02

      Sp2[2, b] += 1.0
    else
      Sp2[3, b] += 1.0
    end

    if splitsupport.C[1, b] == 0x01
      @test splitsupport.C[2, b] == 0x01
      @test splitsupport.C[3, b] == 0x01
      @test splitsupport.C[4, b] == 0x01

      Sp3[1, b] += 1.0
    elseif splitsupport.C[1, b] == 0x02
      @test splitsupport.C[2, b] == 0x02
      @test splitsupport.C[3, b] == 0x02
      @test splitsupport.C[4, b] == 0x02

      Sp3[2, b] += 1.0
    else
      Sp3[3, b] += 1.0
    end

    if brwsupport_1.C[1, b] == 0x01
      @test brwsupport_1.C[2, b] == 0x01

      Sp4[1, b] += 1.0
    elseif brwsupport_1.C[1, b] == 0x02
      @test brwsupport_1.C[2, b] == 0x02

      Sp4[2, b] += 1.0
    else
      Sp4[3, b] += 1.0
    end

    if brwsupport_2.C[1, b] == 0x01
      @test brwsupport_2.C[2, b] == 0x01
      @test brwsupport_2.C[3, b] == 0x01

      Sp5[1, b] += 1.0
    elseif brwsupport_2.C[1, b] == 0x02
      @test brwsupport_2.C[2, b] == 0x02
      @test brwsupport_2.C[3, b] == 0x02

      Sp5[2, b] += 1.0
    else
      Sp5[3, b] += 1.0
    end

    if brwsupport_3.C[1, b] == 0x01
      @test brwsupport_3.C[2, b] == 0x01
      @test brwsupport_3.C[3, b] == 0x01
      @test brwsupport_3.C[4, b] == 0x01

      Sp6[1, b] += 1.0
    elseif brwsupport_3.C[1, b] == 0x02
      @test brwsupport_3.C[2, b] == 0x02
      @test brwsupport_3.C[3, b] == 0x02
      @test brwsupport_3.C[4, b] == 0x02

      Sp6[2, b] += 1.0
    else
      Sp6[3, b] += 1.0
    end
  end

  if all(no.C[no.cl[1:k], :] .== Ctest1)
    Cp1 += 1.0
  end

  if all(mergesupport.C[1:2, :] .== Ctest2)
    Cp2 += 1.0
  end

  if all(splitsupport.C[1:4, :] .== Ctest3)
    Cp3 += 1.0
  end

  if all(brwsupport_1.C[1:2, :] .== Ctest4)
    Cp4 += 1.0
  end

  if all(brwsupport_2.C[1:3, :] .== Ctest5)
    Cp5 += 1.0
  end

  if all(brwsupport_3.C[1:4, :] .== Ctest6)
    Cp6 += 1.0
  end
end

Sp1 /= N
Sp2 /= N
Sp3 /= N
Sp4 /= N
Sp5 /= N
Sp6 /= N

@test maximum(abs(Sp1 - trueSp1)) < 0.005
@test maximum(abs(Sp2 - trueSp2)) < 0.005
@test maximum(abs(Sp3 - trueSp3)) < 0.005
@test maximum(abs(Sp4 - trueSp4)) < 0.005
@test maximum(abs(Sp5 - trueSp5)) < 0.005
@test maximum(abs(Sp6 - trueSp6)) < 0.005

Cp1 /= N
Cp2 /= N
Cp3 /= N
Cp4 /= N
Cp5 /= N
Cp6 /= N

@test_approx_eq_eps Cp1 trueCp1 0.005
@test_approx_eq_eps Cp2 trueCp2 0.005
@test_approx_eq_eps Cp3 trueCp3 0.005
@test_approx_eq_eps Cp4 trueCp4 0.005
@test_approx_eq_eps Cp5 trueCp5 0.005
@test_approx_eq_eps Cp6 trueCp6 0.005

C = copy(no.C)
cl = copy(no.cl)
k = copy(no.k)
v = copy(no.v)
n1s = copy(no.n1s)

logpC = rpostpartitioncols!(C, cl, k, v, n1s, priorC)

@test_approx_eq_eps logpC[1] logpriorC(C, cl, k, priorC) ε
@test_approx_eq_eps logpC[2] logcondpostC(C, cl, k, v, n1s, priorC) ε

cl = [1; 2]
k = 2
v = [2; 4]
n1s = Float64[0 2 0 0; 4 0 3 3]

simcmerge!(k, me.cl[3], me.cl[1], mergevi, mergeni, priorC, mergesupport, me)

@test_approx_eq_eps mergesupport.logpC[1] logpriorC(mergesupport.C, cl, k,
                                                    priorC) ε
@test_approx_eq_eps mergesupport.logpC[2] logcondpostC(mergesupport.C, cl, k, v,
                                                       n1s, priorC) ε

cl = [1; 2; 3; 4]
k = 4
v = [2; 1; 2; 1]
n1s = Float64[0 2 0 0; 1 0 0 1; 2 0 2 2; 1 0 1 0]

simcsplit!(k, sp.cl[1], priorC, splitsupport, sp)

@test_approx_eq_eps splitsupport.logpC[1] logpriorC(splitsupport.C, cl, k,
                                                    priorC) ε
@test_approx_eq_eps splitsupport.logpC[2] logcondpostC(splitsupport.C, cl, k, v,
                                                       n1s, priorC) ε

cl = [1; 2]
k = 2
v = [2; 4]
n1s = Float64[0 2 0 0; 4 0 3 3]

simcbrw!(k, 3, 1, priorC, brwsupport_1, brw_1)

@test_approx_eq_eps brwsupport_1.logpC[1] logpriorC(brwsupport_1.C, cl, k,
                                                    priorC) ε
@test_approx_eq_eps brwsupport_1.logpC[2] logcondpostC(brwsupport_1.C, cl, k, v,
                                                       n1s, priorC) ε

cl = [1; 2; 3]
k = 3
v = [2; 2; 2]
n1s = Float64[0 2 0 0; 2 0 2 2; 2 0 1 1]

simcbrw!(k, 1, 3, priorC, brwsupport_2, brw_2)

@test_approx_eq_eps brwsupport_2.logpC[1] logpriorC(brwsupport_2.C, cl, k,
                                                    priorC) ε
@test_approx_eq_eps brwsupport_2.logpC[2] logcondpostC(brwsupport_2.C, cl, k, v,
                                                       n1s, priorC) ε

cl = [1; 2; 3; 4]
k = 4
v = [2; 1; 2; 1]
n1s = Float64[0 2 0 0; 1 0 0 1; 2 0 2 2; 1 0 1 0]

simcbrw!(k, 1, 4, priorC, brwsupport_3, brw_3)

@test_approx_eq_eps brwsupport_3.logpC[1] logpriorC(brwsupport_3.C, cl, k,
                                                    priorC) ε
@test_approx_eq_eps brwsupport_3.logpC[2] logcondpostC(brwsupport_3.C, cl, k, v,
                                                       n1s, priorC) ε
