output_precision(12)
freq=50;
R1=15e3
R2=5e3
C=0.00002
Vp=230
n=16

vON = 12/19
VS=(Vp/n)-2*vON


%Rectifier
t=linspace(0, 0.2, 1000);
w=2*pi*freq;
vS=VS*cos(w*t);
vstore = zeros(1, length(t));

for i=1:length(t)
  if (vS(i) > 0)
    vstore(i) = vS(i);
  else
    vstore(i) = 0;
  endif
endfor

%envelope detetor
t=linspace(0, 20e-3, 200);
b = cos(w*t)
vS=VS*b;
p = atan(1/(w*R1*C))
tOFF = (1/w) * p
u = exp(-((t-tOFF)/(R1*C)))

vexp = VS * cos(w*tOFF) * u;

vstore = zeros(1, length(t));

for i=1:length(t)
  if (vS(i) > 0)
    vstore(i) = vS(i);
  else
    vstore(i) = 0;
  endif
endfor

f1=figure()
plot(t*1000, vstore)
hold

for i=1:length(t)
  if t(i) < tOFF(1)
    vO(i) = vstore(i);
  elseif vexp(i) > vstore(i)
    vO(i) = vexp(i);
  else 
    vO(i) = vstore(i);
  endif
endfor

plot(t*1000, vO,"y")
title("Full wave Rectifier and Envelope Detector Voltages")
xlabel ("t[ms]")
ylabel ("V[volt]")
legend("rectifier","envelope detector")
print (f1,"rectifier.eps", "-depsc");

f2=figure() 
plot(t*1000, vO,"g")
title("Detailed Envelope Detector")
xlabel ("t[ms]")
ylabel ("V[volt]")
legend("envelope detector")
print (f2,"envelope.eps", "-depsc");

Vavg = (max(vO)+min(vO))/2


%Voltage regulator

vON = 12/19;
R2 = 5000;
ripple = vO - Vavg;

##Ngspice
Is = 5.542088e-4;
Rd = vON/(Is*exp(Vavg/vON));
Rdn = 19*Rd

vOd = (Rdn/(Rdn + R2))*ripple;

Vdiode = 19*vON + vOd;
Vdiode_avg = (max(Vdiode)+min(Vdiode))/2;
Vdeviation = Vdiode_avg - 12

figure
plot(t*1000, Vdiode)
hold
plot(t*1000, Vdiode_avg)
title("Diodes voltage")
xlabel ("t[ms]")
ylabel ("V[volt]")
legend("voltage regulator","average voltage")
print ("diodesaverage.eps", "-depsc");

figure 
plot(t*1000, Vdeviation)
title("Voltage considering deviation")
xlabel ("t[ms]")
ylabel ("V[volt]")
legend("voltage deviation")
print ("deviation.eps", "-depsc");


ripple2 = max(Vdiode) - min(Vdiode);
figure 
plot(t*1000, ripple2)
title("Voltage Ripple")
xlabel ("t[ms]")
ylabel ("V[volt]")
legend("Ripple")
print ("ripple.eps", "-depsc");





%figure of merit


Cost_R = 20;
Cost_C = 20;
Cost_Dio = 0.1 * 23;
cost = Cost_R + Cost_C + Cost_Dio;
Merit = 1 ./ ( cost * (ripple2 + Vdeviation + 10e-6));


fp = fopen ("merit.tex", "w");

fprintf(fp, "Total Cost & %e \\\\ \\hline \n", cost);
fprintf(fp, "Merit & %e \\\\ \n", Merit);
fclose (fp);

fp2 = fopen ("ripple.tex", "w");
fprintf(fp2, "Ripple & %e \\\\ \\hline \n", ripple2);
fclose (fp2);

