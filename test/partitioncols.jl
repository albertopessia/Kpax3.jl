# This file is part of Kpax3. License is MIT.

ε = 1.0e-14

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
  ω = [1.0; 1.0; 1.0 - 1.0 / k; 1.0 / k]

  for γ in ([1.0; 0.0; 0.0], [0.0; 1.0; 0.0], [0.0; 0.0; 1.0],
            [0.4; 0.3; 0.3], [0.5; 0.3; 0.2], [0.7; 0.2; 0.1],
            [0.1; 0.1; 0.1], [0.3; 0.1; 0.1], [0.0; 0.2; 0.1])
    x1 = AminoAcidPriorCol(data, k, γ, r1)
    x2 = AminoAcidPriorCol(data, k, γ, r2)

    for s in 1:3
      @test_approx_eq_eps x1.logγ[s] log(γ[s] / sum(γ)) ε
      @test_approx_eq_eps x1.logω[s] log(ω[s]) ε

      @test_approx_eq_eps x2.logγ[s] log(γ[s] / sum(γ)) ε
      @test_approx_eq_eps x2.logω[s] log(ω[s]) ε
    end

    for b in 1:m, s in 1:4
      @test_approx_eq_eps x1.A[s, b] A1[s, b] ε
      @test_approx_eq_eps x1.B[s, b] B1[s, b] ε

      @test_approx_eq_eps x2.A[s, b] A2[s, b] ε
      @test_approx_eq_eps x2.B[s, b] B2[s, b] ε
    end
  end
end

settings = KSettings("typesmcmc.bin", 1, 0, 1, [1.0; 0.0; 0.0], 0.0, 1.0,
                     [0.6; 0.35; 0.05], 135.0, 1.0, 1.0, 5.0, 100, 1, true, 1)

data = UInt8[0x01 0x01 0x01 0x00 0x00 0x00;
             0x00 0x00 0x00 0x01 0x01 0x00;
             0x01 0x01 0x00 0x01 0x00 0x01;
             0x01 0x01 0x00 0x00 0x00 0x00]

m, n = size(data)

R = [13; 13; 13; 42; 42; 76]
k = length(unique(R))

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)

no = AminoAcidMCMC(data, R, priorR, priorC, settings)

C = zeros(UInt8, settings.maxclust, m)

cs = ([0x01; 0x01; 0x01], [0x02; 0x02; 0x02], [0x03; 0x03; 0x03],
      [0x03; 0x03; 0x04], [0x03; 0x04; 0x03], [0x04; 0x03; 0x03],
      [0x03; 0x04; 0x04], [0x04; 0x03; 0x04], [0x04; 0x04; 0x03],
      [0x04; 0x04; 0x04])

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

  logp1[l] = logpriorC(no.C, no.cl, priorC.logγ, priorC.logω)
  logp2[l] = logpriorC(C, k, priorC.logγ, priorC.logω)
  logp3[l] = logcondpostC(no.C, no.cl, no.v, no.n1s, priorC.logω, priorC)
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

ss = [0x01; 0x02; 0x03]
logp4 = zeros(Float64, 3^m)
l = 0

for s1 in ss, s2 in ss, s3 in ss, s4 in ss
  l += 1
  S = [s1; s2; s3; s4]
  logp4[l] = logcondpostS(S, no.cl, no.v, no.n1s, priorC.logω, priorC)
end

M = maximum(logp4)
p4 = exp(M + log(sum(exp(logp4 - M))))

@test_approx_eq_eps p4 1.0 ε

# suppose we merge cluster 2 and cluster 3
me = AminoAcidMCMC(data, R, priorR, priorC, settings)

mergesupport = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))
mergevi = me.v[me.cl[2]] + me.v[me.cl[3]]
mergeni = vec(no.n1s[no.cl[2], :] + no.n1s[no.cl[3], :])
mergesupport.logω = [0.0; 0.0; log(k - 2.0) - log(k - 1.0); - log(k - 1.0)]

# suppose we split cluster 1
sp = AminoAcidMCMC(data, R, priorR, priorC, settings)

splitsupport = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))
splitsupport.vi = sp.v[sp.cl[1]] - 1
splitsupport.ni = vec(sp.n1s[sp.cl[1], :]) - float(data[:, 3])
splitsupport.vj = 1
splitsupport.nj = float(data[:, 3])
splitsupport.logω = [0.0; 0.0; log(k) - log(k + 1.0); - log(k + 1.0)]

# suppose we move unit 6 into cluster 2 (merge)
brw_1 = AminoAcidMCMC(data, R, priorR, priorC, settings)
brwsupport_1 = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))

brwsupport_1.ni = float(data[:, 6])
brwsupport_1.logω = [0.0; 0.0; log(k - 2.0) - log(k - 1.0); - log(k - 1.0)]

# suppose we move unit 6 into the same cluster
brw_2 = AminoAcidMCMC(data, R, priorR, priorC, settings)
brwsupport_2 = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))

brwsupport_2.ni = float(data[:, 6])
brwsupport_2.logω = [0.0; 0.0; log(k - 1.0) - log(k); - log(k)]

# suppose we move unit 3 into cluster 2
brw_3 = AminoAcidMCMC(data, R, priorR, priorC, settings)
brwsupport_3 = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))

brwsupport_3.ni = float(data[:, 3])
brwsupport_3.logω = [0.0; 0.0; log(k - 1.0) - log(k); - log(k)]

# suppose we move unit 3 into its own cluster (split)
brw_4 = AminoAcidMCMC(data, R, priorR, priorC, settings)
brwsupport_4 = KSupport(size(data, 1), size(data, 2), settings.maxclust,
                        size(data, 2))

brwsupport_4.ni = float(data[:, 3])
brwsupport_4.logω = [0.0; 0.0; log(k) - log(k + 1.0); - log(k + 1.0)]

N = 10000000

Sp1 = zeros(Float64, 3, m)
Sp2 = zeros(Float64, 3, m)
Sp3 = zeros(Float64, 3, m)
Sp4 = zeros(Float64, 3, m)
Sp5 = zeros(Float64, 3, m)
Sp6 = zeros(Float64, 3, m)
Sp7 = zeros(Float64, 3, m)

Cp1 = zero(Float64)
Cp2 = zero(Float64)
Cp3 = zero(Float64)
Cp4 = zero(Float64)
Cp5 = zero(Float64)
Cp6 = zero(Float64)
Cp7 = zero(Float64)

Ctest1 = [0x04 0x02 0x01 0x02; 0x03 0x02 0x01 0x02; 0x03 0x02 0x01 0x02]
Ctest2 = [0x04 0x02 0x01 0x02; 0x03 0x02 0x01 0x02]
Ctest3 = [0x02 0x03 0x01 0x02; 0x02 0x04 0x01 0x02; 0x02 0x03 0x01 0x02;
          0x02 0x03 0x01 0x02]
Ctest4 = [0x04 0x01 0x01 0x02; 0x03 0x01 0x01 0x02]
Ctest5 = [0x01 0x03 0x02 0x01; 0x01 0x04 0x02 0x01; 0x01 0x03 0x02 0x01]
Ctest6 = [0x01 0x02 0x01 0x04; 0x01 0x02 0x01 0x03; 0x01 0x02 0x01 0x03]
Ctest7 = [0x01 0x01 0x02 0x04; 0x01 0x01 0x02 0x03; 0x01 0x01 0x02 0x03;
          0x01 0x01 0x02 0x03]

#=
R = [13; 13; 80; 42; 42; 76]
k = length(unique(R))

priorR = EwensPitman(settings.α, settings.θ)
priorC = AminoAcidPriorCol(data, k, settings.γ, settings.r)

ss = [0x01; 0x02; 0x03]
obj = AminoAcidMCMC(data, R, priorR, priorC, settings)

trueSp = zeros(Float64, 3, 4)
for s1 in ss, s2 in ss, s3 in ss, s4 in ss
  S = [s1; s2; s3; s4]
  p = exp(logcondpostS(S, obj.cl, obj.v, obj.n1s, priorC.logω, priorC))

  trueSp[s1, 1] += p
  trueSp[s2, 2] += p
  trueSp[s3, 3] += p
  trueSp[s4, 4] += p
end

for i in 1:12
  @printf("%.100f\n", trueSp[i])
end

obj.C[obj.cl, :] = Ctest7

p = exp(logcondpostC(obj.C, obj.cl, obj.v, obj.n1s, priorC.logω, priorC))
@printf("%.100f\n", p)
=#

trueSp1 = hcat([0.400605930740163984626889259743620641529560089111328125000000
                0.402853999936236828460778269800357520580291748046875000000000
                0.196540069323599270179059317342762369662523269653320312500000],
               [0.469929717191958196131906788650667294859886169433593750000000
                0.356242671218090556362056986472452990710735321044921875000000
                0.173827611589951080972582531103398650884628295898437500000000],
               [0.796475635596692543849428602698026224970817565917968750000000
                0.203499409017467020044378500642778817564249038696289062500000
                0.000024955385840355038363580844618105913923500338569283485413],
               [0.730452799237661709597091430623549968004226684570312500000000
                0.266615119904645703208245777204865589737892150878906250000000
                0.002932080857692189075625055494356274721212685108184814453125])

trueSp2 = hcat([0.343752696310923200329057181079406291246414184570312500000000
                0.424245751920919844657476005522767081856727600097656250000000
                0.232001551768157093791344891542394179850816726684570312500000],
               [0.669677308216207967106470277940388768911361694335937500000000
                0.326357115550276744020408159485668875277042388916015625000000
                0.003965576233515530867046461338532026275061070919036865234375],
               [0.788118240486919030551860032574040815234184265136718750000000
                0.211824568291017867327497015139670111238956451416015625000000
                0.000057191222062770490303443282620321497233817353844642639160],
               [0.669677308216207967106470277940388768911361694335937500000000
                0.326357115550276744020408159485668875277042388916015625000000
                0.003965576233515523928152557431303648627363145351409912109375])

trueSp3 = hcat([0.508985055567922284325277360039763152599334716796875000000000
                0.417055884580573688058535708478302694857120513916015625000000
                0.073959059851504208027428433069871971383690834045410156250000],
               [0.536479879971808792937792986776912584900856018066406250000000
                0.304601180625598655371533141078543849289417266845703125000000
                0.158918939402592551690673872144543565809726715087890625000000],
               [0.714885118733651014899521669576643034815788269042968750000000
                0.284126957067931396050397552244248799979686737060546875000000
                0.000987924198417838416580449845127986918669193983078002929688],
               [0.536479879971808792937792986776912584900856018066406250000000
                0.304601180625598655371533141078543849289417266845703125000000
                0.158918939402592440668371409628889523446559906005859375000000])

trueSp4 = hcat([0.343752696310923089306754718563752248883247375488281250000000
                0.424245751920919733635173543007113039493560791015625000000000
                0.232001551768157121546920507171307690441608428955078125000000],
               [0.669677308216208189151075202971696853637695312500000000000000
                0.326357115550276744020408159485668875277042388916015625000000
                0.003965576233515533469131675303742667892947793006896972656250],
               [0.788118240486919474641069882636656984686851501464843750000000
                0.211824568291018006105375093284237664192914962768554687500000
                0.000057191222062770483527179704585918784687237348407506942749],
               [0.669677308216208300173377665487350896000862121582031250000000
                0.326357115550276688509256928227841854095458984375000000000000
                0.003965576233515524795514295419707195833325386047363281250000])

trueSp5 = hcat([0.400605930740163707071133103454485535621643066406250000000000
                0.402853999936236606416173344769049435853958129882812500000000
                0.196540069323599075890030007940367795526981353759765625000000],
               [0.469929717191957918576150632361532188951969146728515625000000
                0.356242671218090334317452061441144905984401702880859375000000
                0.173827611589950969950280068587744608521461486816406250000000],
               [0.796475635596692099760218752635410055518150329589843750000000
                0.203499409017466853510924806869297754019498825073242187500000
                0.000024955385840355021422921899532099132557050324976444244385],
               [0.730452799237661598574788968107895925641059875488281250000000
                0.266615119904645647697094545947038568556308746337890625000000
                0.002932080857692187774582448511750953912269324064254760742188])

trueSp6 = hcat([0.652563188862539500512127688125474378466606140136718750000000
                0.343736988372037033379058357240865007042884826660156250000000
                0.003699822765423427944897483143904537428170442581176757812500],
               [0.730452799237662042663998818170512095093727111816406250000000
                0.266615119904645703208245777204865589737892150878906250000000
                0.002932080857692189509305924488558048324193805456161499023438],
               [0.731536451339322657538843941438244655728340148925781250000000
                0.267010652696526185057024349589482881128787994384765625000000
                0.001452895964151075221607034571036365377949550747871398925781],
               [0.469929717191958362665360482424148358404636383056640625000000
                0.356242671218090389828603292698971927165985107421875000000000
                0.173827611589951247506036224876879714429378509521484375000000])

trueSp7 = hcat([0.508985055567922173302974897524109110236167907714843750000000
                0.417055884580573688058535708478302694857120513916015625000000
                0.073959059851504208027428433069871971383690834045410156250000],
               [0.536479879971808903960095449292566627264022827148437500000000
                0.304601180625598766393835603594197891652584075927734375000000
                0.158918939402592634957400719031284097582101821899414062500000],
               [0.714885118733650792854916744545334950089454650878906250000000
                0.284126957067931451561548783502075821161270141601562500000000
                0.000987924198417839067101753336430647323140874505043029785156],
               [0.536479879971809014982397911808220669627189636230468750000000
                0.304601180625598877416138066109851934015750885009765625000000
                0.158918939402592523935098256515630055218935012817382812500000])

trueCp1 = 0.0148123191595730674396946824344922788441181182861328125
trueCp2 = 0.0194745023689023473434378530555477482266724109649658203125
trueCp3 = 0.01435615231057929368219117094440662185661494731903076171875
trueCp4 = 0.0399612317423107266112225488541298545897006988525390625
trueCp5 = 0.0103107627899840618990179308411825331859290599822998046875
trueCp6 = 0.022037454503308018249896349516347981989383697509765625
trueCp7 = 0.01226441605595194640765388527370305382646620273590087890625

for t in 1:N
  rpostpartitioncols!(no.C, no.cl, no.v, no.n1s, priorC)

  simcmerge!(k - 1, me.cl[2], me.cl[3], mergevi, mergeni, priorC,
             mergesupport, me)

  simcsplit!(k + 1, sp.cl[1], priorC, splitsupport, sp)

  simcbrw!(k - 1, 76, 42, priorC, brwsupport_1, brw_1)
  simcbrw!(k, 76, 76, priorC, brwsupport_2, brw_2)
  simcbrw!(k, 13, 42, priorC, brwsupport_3, brw_3)
  simcbrw!(k + 1, 13, 80, priorC, brwsupport_4, brw_4)

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

      Sp6[1, b] += 1.0
    elseif brwsupport_3.C[1, b] == 0x02
      @test brwsupport_3.C[2, b] == 0x02
      @test brwsupport_3.C[3, b] == 0x02

      Sp6[2, b] += 1.0
    else
      Sp6[3, b] += 1.0
    end

    if brwsupport_4.C[1, b] == 0x01
      @test brwsupport_4.C[2, b] == 0x01
      @test brwsupport_4.C[3, b] == 0x01
      @test brwsupport_4.C[4, b] == 0x01

      Sp7[1, b] += 1.0
    elseif brwsupport_4.C[1, b] == 0x02
      @test brwsupport_4.C[2, b] == 0x02
      @test brwsupport_4.C[3, b] == 0x02
      @test brwsupport_4.C[4, b] == 0x02

      Sp7[2, b] += 1.0
    else
      Sp7[3, b] += 1.0
    end
  end

  if all(no.C[no.cl, :] .== Ctest1)
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

  if all(brwsupport_3.C[1:3, :] .== Ctest6)
    Cp6 += 1.0
  end

  if all(brwsupport_4.C[1:4, :] .== Ctest7)
    Cp7 += 1.0
  end
end

Sp1 /= N
Sp2 /= N
Sp3 /= N
Sp4 /= N
Sp5 /= N
Sp6 /= N
Sp7 /= N

@test maximum(abs(Sp1 - trueSp1)) < 0.0005
@test maximum(abs(Sp2 - trueSp2)) < 0.0005
@test maximum(abs(Sp3 - trueSp3)) < 0.0005
@test maximum(abs(Sp4 - trueSp4)) < 0.0005
@test maximum(abs(Sp5 - trueSp5)) < 0.0005
@test maximum(abs(Sp6 - trueSp6)) < 0.0005
@test maximum(abs(Sp7 - trueSp7)) < 0.0005

Cp1 /= N
Cp2 /= N
Cp3 /= N
Cp4 /= N
Cp5 /= N
Cp6 /= N
Cp7 /= N

@test_approx_eq_eps Cp1 trueCp1 0.0005
@test_approx_eq_eps Cp2 trueCp2 0.0005
@test_approx_eq_eps Cp3 trueCp3 0.0005
@test_approx_eq_eps Cp4 trueCp4 0.0005
@test_approx_eq_eps Cp5 trueCp5 0.0005
@test_approx_eq_eps Cp6 trueCp6 0.0005
@test_approx_eq_eps Cp7 trueCp7 0.0005

C = copy(no.C)
cl = copy(no.cl)
v = copy(no.v)
n1s = copy(no.n1s)
logω = copy(priorC.logω)

logpC = rpostpartitioncols!(C, cl, v, n1s, priorC)

@test_approx_eq_eps logpC[1] logpriorC(C, cl, priorC.logγ, logω) ε
@test_approx_eq_eps logpC[2] logcondpostC(C, cl, v, n1s, logω, priorC) ε

cl = [1; 2]
v = [3; 3]
n1s = [3.0 0.0 2.0 2.0; 0.0 2.0 2.0 0.0]
logω = copy(mergesupport.logω)

simcmerge!(k - 1, me.cl[2], me.cl[3], mergevi, mergeni, priorC, mergesupport,
           me)

@test_approx_eq_eps mergesupport.logpC[1] logpriorC(mergesupport.C, cl,
                                                    priorC.logγ, logω) ε
@test_approx_eq_eps mergesupport.logpC[2] logcondpostC(mergesupport.C, cl, v,
                                                       n1s, logω, priorC) ε

cl = [1; 2; 3; 4]
v = [2; 2; 1; 1]
n1s = [2.0 0.0 2.0 2.0; 0.0 2.0 1.0 0.0; 0.0 0.0 1.0 0.0; 1.0 0.0 0.0 0.0]
logω = copy(splitsupport.logω)

simcsplit!(k + 1, sp.cl[1], priorC, splitsupport, sp)

@test_approx_eq_eps splitsupport.logpC[1] logpriorC(splitsupport.C, cl,
                                                    priorC.logγ, logω) ε
@test_approx_eq_eps splitsupport.logpC[2] logcondpostC(splitsupport.C, cl, v,
                                                       n1s, logω, priorC) ε

cl = [1; 2]
v = [3; 3]
n1s = [3.0 0.0 2.0 2.0; 0.0 2.0 2.0 0.0]
logω = copy(brwsupport_1.logω)

simcbrw!(k - 1, 76, 42, priorC, brwsupport_1, brw_1)

@test_approx_eq_eps brwsupport_1.logpC[1] logpriorC(brwsupport_1.C, cl,
                                                    priorC.logγ, logω) ε
@test_approx_eq_eps brwsupport_1.logpC[2] logcondpostC(brwsupport_1.C, cl, v,
                                                       n1s, logω, priorC) ε

cl = [1; 2; 3]
v = [3; 2; 1]
n1s = [3.0 0.0 2.0 2.0; 0.0 2.0 1.0 0.0; 0.0 0.0 1.0 0.0]
logω = copy(brwsupport_2.logω)

simcbrw!(k, 76, 76, priorC, brwsupport_2, brw_2)

@test_approx_eq_eps brwsupport_2.logpC[1] logpriorC(brwsupport_2.C, cl,
                                                    priorC.logγ, logω) ε
@test_approx_eq_eps brwsupport_2.logpC[2] logcondpostC(brwsupport_2.C, cl, v,
                                                       n1s, logω, priorC) ε

cl = [1; 2; 3]
v = [2; 3; 1]
n1s = [2.0 0.0 2.0 2.0; 1.0 2.0 1.0 0.0; 0.0 0.0 1.0 0.0]
logω = copy(brwsupport_3.logω)

simcbrw!(k, 13, 42, priorC, brwsupport_3, brw_3)

@test_approx_eq_eps brwsupport_3.logpC[1] logpriorC(brwsupport_3.C, cl,
                                                    priorC.logγ, logω) ε
@test_approx_eq_eps brwsupport_3.logpC[2] logcondpostC(brwsupport_3.C, cl, v,
                                                       n1s, logω, priorC) ε

cl = [1; 2; 3; 4]
v = [2; 2; 1; 1]
n1s = [2.0 0.0 2.0 2.0; 0.0 2.0 1.0 0.0; 0.0 0.0 1.0 0.0; 1.0 0.0 0.0 0.0]
logω = copy(brwsupport_4.logω)

simcbrw!(k + 1, 13, 80, priorC, brwsupport_4, brw_4)

@test_approx_eq_eps brwsupport_4.logpC[1] logpriorC(brwsupport_4.C, cl,
                                                    priorC.logγ, logω) ε
@test_approx_eq_eps brwsupport_4.logpC[2] logcondpostC(brwsupport_4.C, cl, v,
                                                       n1s, logω, priorC) ε
