R1=1000;
R2=1000;
R3=400000;
R4=1000;
C1=0.0733e-6;
C2=0.22e-6;

freq=logspace(1,8,100);

w=2*pi*freq;
S = w*j;



% Ts = (R1*C1*2*pi*freq*j)/(1+R1*C1*2*pi*freq*j)(1 + (R3/R4))(1/(1+R2*C2*2*pi*freq*j))

wL = 1/(R2*C2)
wH = 1/(R1*C1)
wO = sqrt(wL*wH)
fO=wO/(2*pi)

gain =abs(((R1*C1*wO*j)/(1+R1*C1*wO*j))*(1+(R3/R4))*(1/(1+R2*C2*wO*j)))
gaindb = 20*log10(abs(gain))

z_in = (R1 + 1/(j*wO*C1))
z_out = R2/(j*wO*C2)/ (R2+1/(j*wO*C2))

cost =13726.88
gain_dev= abs(gain-100)
central_freq_deviation =  abs(fO-1000) 
Merit = 1/(cost*(gain_dev+central_freq_deviation +1e-6))

fp=fopen('../doc/valintro.txt',"w");
fprintf(fp,"|R1| & %f \\\\ \\hline\n|R2| & %f \\\\ \\hline\n|R3| & %f \\\\ \\hline\n|R4| & %f \\\\ \\hline\n|C1| & %f \\\\ \\hline\n|C2| & %f \\\\ \\hline",R1,R2,R3,R4,C1,C2);
fclose(fp);

fp=fopen('../doc/a1.txt',"w");
fprintf(fp,"|Central Frequency| & %f \\\\ \\hline\n|Gain| & %f \\\\ \\hline",fO,gain);
fclose(fp);

fp=fopen('../doc/a2.txt',"w");
fprintf(fp,"|Input Impedance| & %f \\\\ \\hline\n|Output Impedance| & %f \\\\ \\hline",z_in,z_out);

fp=fopen('../doc/a3.txt',"w");
fprintf(fp,"|Cost| & %f \\\\ \\hline\n|Merit| & %f \\\\ \\hline",cost,Merit);
fclose(fp);

Ts = (R1*C1*2*pi*freq*j)./(1+R1*C1*2*pi*freq*j)*(1+R3/R4).*(1./(1+R2*C2*2*pi*freq*j));

f1 = figure();
semilogx(freq,20*log10(abs(Ts)));
xlabel('Frequency [Hz]');
ylabel('Gain[dB]')
print(f1, 'gain.eps', '-depsc')

f2 = figure();
semilogx(freq,180*arg(Ts)/pi);
xlabel('Frequency[Hz]');
ylabel('Phase[rad]')
print(f2, 'phase.eps', '-depsc')


