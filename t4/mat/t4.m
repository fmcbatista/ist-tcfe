VT=25e-3
BFN=178.7
VAFN=69.7
RE=200
RC=1100
R1=55000
R2=10000
VBEON=0.7
VCC=12
Rin=100
Ci=0.1e-3
Cb=2.8e-3




Req=1/(1/R1+1/R2)
VEQ=R2/(R1+R2)*VCC
IB1=(VEQ-VBEON)/(Req+(1+BFN)*RE)
IC1=BFN*IB1
IE1=(1+BFN)*IB1
VE1=RE*IE1
VO1=VCC-RC*IC1
VCE=VO1-VE1


gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

RSB=Req*Rin/(Req+Rin)

AV1 = RSB/Rin * RC*(RE-gm1*rpi1*ro1)/((ro1+RC+RE)*(RSB+rpi1+RE)+gm1*RE*ro1*rpi1 - RE^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple = Req/(Req+Rin) * gm1*RC/(1+gm1*RE)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE=0

ZI1s= 1/((1/Req)+(1/rpi1))
ZO1s= 1/(1/ro1+1/RC)

AV1 = gm1*ZO1s*(ZI1s/(ZI1s+Rin))
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/Rin * gm1*RC/(1+gm1*RE)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE=73.8
ZI1= 1/(1/Req+1/(((ro1+RC+RE)*(rpi1+RE)+gm1*RE*ro1*rpi1 - RE^2)/(ro1+RC+RE)))
ZX = ro1*((RSB+rpi1)*RE/(RSB+rpi1+RE))/(1/(1/ro1+1/(rpi1+RSB)+1/RE+gm1*rpi1/(rpi1+RSB)))
ZX = ro1*(   1/RE+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC)


%ouput stage
BFP = 227.3
VAFP = 37.2
Rout = 220
VEBON = 0.7
Co=1.52e-3
RL=8

VI2 = VO1
IE2 = (VCC-VEBON-VI2)/Rout
IC2 = BFP/(BFP+1)*IE2
IB2=IE2-IC2
VO2 = VCC - Rout*IE2

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/Rout

AV2 = gm2/(gm2+gpi2+go2+ge2)
AVDB2 = 20*log10(abs(AV2))
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)


%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

freq=logspace(1,8,100);

w = 2*pi*freq;

G=1:length(freq)

for k=1:length(freq)

ZCi(k) = 1/(j*w(k)*Ci);
ZCb(k) = 1/(j*w(k)*Cb);
ZCo(k) = 1/(j*w(k)*Co);

Zin(k) = ZCi(k) + Rin;
ZB(k) = 1/((1/R1)+(1/R2));
Zbyp(k) = 1/((1/RE)+(1/ZCb(k)));

gpi1=1/rpi1
Yin(k) = 1/Zin(k);
YB(k) = 1/ZB(k);
Ybyp(k) = 1/Zbyp(k);
go1 = 1/ro1
gc = 1/RC
gout = 1/Rout
gE = 1/RE 

A = [0 , 0 , -gpi2-gm2 , gout+gpi2+gm2+go2+(1/(ZCo(k))) , -1/(ZCo(k)) , 0 ; gm1 , -gm1 , go1+gc+gpi2 , -go1-gpi2 , 0 , 0 ; -gpi1-gm1 , gE+gpi1+gm1+go1+(1/(ZCb(k))) , -go1 , 0 , 0 , 0 ; (1/(ZCi(k)))+gpi1 , -gpi1 , 0 , 0 , 0 , -1/(ZCi(k)) ; -1/(ZCi(k)) , 0 , 0 , 0 , 0 , (1/Rin)+(1/(ZCi(k))) ; 0 , 0 , 0 , -1/(ZCo(k)) , (1/(ZCo(k))) + (1/RL) , 0];

B = [0 ; 0 ; 0 ; 0 ; 0.01/Rin ; 0];

X = inv(A);

S = (inv(A)) * B;

G(k)=S(5)/S(6);

endfor

Gdb = 20*log10(abs(G));

%Gain = max(abs(G))
%Gaindb = max(Gdb)
%y1=(max(Gdb))-3
%LowerCutoffFrequency=interp1(Gdb,freq,y1)

f1=figure()
plot(log10(freq),Gdb)
xlabel("Frequency [Hz]")
ylabel("Gain [dB]")
title("Frequency Response")
print(f1,"freqresponse.eps","-depsc")


fp = fopen('../doc/TA1.tex',"w");
fprintf(fp,"|VT| & %f \\\\ \\hline\n|beta| & %f \\\\ \\hline\n|VA| & %f \\\\ \\hline\n|VBEON| & %f \\\\ \\hline\n|Vcc| & %f \\\\ \\hline\n|Req| & %f \\\\ \\hline\n|Veq| & %f \\\\ \\hline\n|IB1| & %f \\\\ \\hline\n|IC1| & %f \\\\ \\hline\n|IE1| & %f \\\\ \\hline\n |$V_{emit}$| & %f \\\\ \\hline\n|$V_{coll}$| & %f  \\\\ \\hline\n|VCE| & %f \\\\ \\hline|$V_{base}$| & %f  \\\\ \\hline" , VT,BFN,VAFN,VBEON,VCC,Req,VEQ,IB1,IC1,IE1,VE1,VO1,VCE,(VBEON+VE1));
fclose(fp);

fp = fopen('../doc/TA2.tex',"w");
fprintf(fp,"|gm1| & %f \\\\ \\hline\n|ro1| & %f \\\\ \\hline\n|rpi| & %f \\\\ \\hline",gm1,ro1,rpi1);
fclose(fp);


fp = fopen('../doc/TA3.tex',"w");
fprintf(fp,"|Input impedance| & %f \\\\ \\hline\n|Output impedance| & %f \\\\ \\hline\n|Gain| & %f \\\\ \\hline\n|Gain(dB)| & %f \\\\ \\hline",ZI1s,ZO1s,AV1,AVI_DB);
fclose(fp);


fp = fopen('../doc/TA4.tex',"w");
fprintf(fp,"|VO2| & %f \\\\ \\hline\n|$V_{coll}$| & %f \\\\ \\hline\n|VA| & %f \\\\ \\hline\n|$V_{emit}$| & %f \\\\ \\hline\n|Vcc| & %f \\\\ \\hline\n|IB2| & %f \\\\ \\hline\n|IC2| & %f \\\\ \\hline\n|IE2| & %f \\\\ \\hline",VO2,VO1,VO2,VCC,IB2,IC2,IE2);
fclose(fp);

fp = fopen('../doc/TA5.tex',"w");
fprintf(fp,"|gm2| & %f \\\\ \\hline\n|go2| & %f \\\\ \\hline\n|gpi2| & %f \\\\ \\hline\n|ge2| & %f \\\\ \\hline\n|Input Impedance| & %f \\\\ \\hline\n|Output Impedance| & %f \\\\ \\hline\n|Gain| & %f \\\\ \\hline\n|Gain(dB)| & %f \\\\ \\hline",gm2,go2,gpi2,ge2,ZI2,ZO2,AV2,AVDB2);
fclose(fp);




