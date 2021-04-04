pkg load control
output_precision(10)
format long
fp = fopen('Values.txt',"r");
values = dlmread(fp,'\n');
R1 = values(1)
R2 = values(2)
R3 = values(3)
R4 = values(4)
R5 = values(5)
R6 = values(6)
R7 = values(7)
Vs = values(8)
Kb = values(9)
Kd = values(10)
C = values(11)
null=0;
A = [1,0,0,0,0,0,0 ; -(1/R1),(1/R3)+(1/R2)+(1/R1),-(1/R2),-(1/R3),0,0,0 ; 0,(Kb+(1/R2)),(-1/R2),(-Kb),0,0,0 ; 0,(1/R3),0,(-(1/R4)-(1/R5)-(1/R3)),(1/R5),(1/R7),-(1/R7) ; 0,(-Kb),0,(Kb+(1/R5)),(-1/R5),0,0 ; 0,0,0,0,0,((1/R6)+(1/R7)),(-1/R7);0,0,0,1,0,Kd/R6,-1];

B = inv(A);

E = [Vs;0;0;0;0;0;0];

D = B*E;

fp=fopen('/home/fmcb/ist-tcfe/t2/doc/TA1.tex',"w");
fprintf(fp,"v1 & %f \\\\ \\hline\nv2 & %f \\\\ \\hline\nv3 & %f \\\\ \\hline\nv4 & %f \\\\ \\hline\nv5 & %f \\\\ \\hline\nv6 & %f \\\\ \\hline\nv7 & %f \\\\ \\hline\nv8 & %f \\\\ \\hline",D(1),D(2),D(3),null,D(4),D(5),D(6),D(7));
fclose(fp);

fp2 = fopen('ngspice_values.txt',"w");
fprintf(fp2, "R1 1 2 %f\n",values(1));
fprintf(fp2, "R2 3 2 %f\n",values(2));
fprintf(fp2, "R3 2 5 %f\n",values(3));
fprintf(fp2, "R4 GND 5 %f\n",values(4));
fprintf(fp2, "R5 5 6 %f\n",values(5));
fprintf(fp2, "R6 GND 4 %f\n",values(6));
fprintf(fp2, "R7 7 8 %f\n",values(7));
fprintf(fp2, "C 6 8 %f\n",values(11));
fprintf(fp2, "Vs 1 GND %f\n",values(8));
fprintf(fp2, "V3 4 7 DC 0 \n");
fprintf(fp2, "Hc 5 8 V3 %f\n",values(10));
fprintf(fp2, "Gb 6 3 (2,5) %f\n",values(9));
fclose(fp2);


%DADOS PARA O NGSPICE 2

fp3 = fopen('ngspice_2.txt',"w");
fprintf(fp3, "R1 GND 2 %f\n",values(1));
fprintf(fp3, "R2 3 2 %f\n",values(2));
fprintf(fp3, "R3 2 5 %f\n",values(3));
fprintf(fp3, "R4 GND 5 %f\n",values(4));
fprintf(fp3, "R5 5 6 %f\n",values(5));
fprintf(fp3, "R6 GND 4 %f\n",values(6));
fprintf(fp3, "R7 7 8 %f\n",values(7));
fprintf(fp3, "Vx 6 8 %f\n",D(5)-D(7));
fprintf(fp3, "Vs 1 GND 0 \n");
fprintf(fp3, "V3 4 7 DC 0 \n");
fprintf(fp3, "Hc 5 8 V3 %f\n",values(10));
fprintf(fp3, "Gb 6 3 (2,5) %f\n",values(9));

fclose(fp3);

%theoretical analysis 2

H = [-(1/R1),0, -(1/R4),0, -(1/R6),0,0,0,0;-(1/R1)-(1/R2)-(1/R3), 1/R2,1/R3,0,0,0,0,0,0; (-Kb -(1/R2)),1/R2,Kb,0,0,0,0,0,0;0,0,1,0,Kd/R6,-1,0,0,0;0,0,0,0,(1/R6+1/R7),-(1/R7),0,0,0;0,0,0,1,0,-1,0,0,0;Kb,0,-Kb,0,0,0,0,-1,0;0,0,1/R5,-(1/R5),0,0,-1,0,0;0,0,0,0,0,0,-1,1,-1];

Vx = D(5)-D(7)

Z= [0;0;0;0;0;Vx;0;0;0];

T = inv(H);

X = T*Z;

Req = (X(4)-X(6))/X(9);
sprintf("%.6f\n",X);

fp=fopen('/home/fmcb/ist-tcfe/t2/doc/TA2.tex',"w");
fprintf(fp,"v1 & %f \\\\ \\hline\nv2 & %f \\\\ \\hline\nv3 & %f \\\\ \\hline\nv4 & %f \\\\ \\hline\nv5 & %f \\\\ \\hline\nv6 & %f \\\\ \\hline\nv7 & %f \\\\ \\hline\nv8 & %f \\\\ \\hline\nIx & %f \\\\ \\hline\nVx & %f \\\\ \\hline\nReq & %f \\\\ \\hline",null,X(1),X(2),null,X(3),X(4),X(5),X(6),,,Req);
fclose(fp);


%theoretical analysis 3
fp4 = fopen('ngspice_3.txt',"w");
fprintf(fp4, "R1 1 2 %f\n",values(1));
fprintf(fp4, "R2 3 2 %f\n",values(2));
fprintf(fp4, "R3 2 5 %f\n",values(3));
fprintf(fp4, "R4 GND 5 %f\n",values(4));
fprintf(fp4, "R5 5 6 %f\n",values(5));
fprintf(fp4, "R6 GND 4 %f\n",values(6));
fprintf(fp4, "R7 7 8 %f\n",values(7));
fprintf(fp4, "C 6 8 %f\n",values(11));
fprintf(fp4, "Vs 1 GND 0 \n");
fprintf(fp4, "V3 4 7 DC 0 \n");
fprintf(fp4, "Hc 5 8 V3 %f\n",values(10));
fprintf(fp4, "Gb 6 3 (2,5) %f\n",values(9));

fclose(fp4);


Tau = Req*C;

t = 0:2e-6:20e-3;

v6n = X(4)*exp(-t/Tau);

v8n = X(6)*exp(-t/Tau);

f1 = figure();

plot(t*1000,v6n,"b");

xlabel("t[ms]");

ylabel("v6n[V]");

title("Natural solution");

print(f1,"natural_solution.eps","-depsc");


%forced solution

w = 2*pi*1000;
Zc = 1/(j*w*C);

G = [1,0,0,0,0,0,0; 1/R1, -(1/R1)-(1/R2)-(1/R3), 1/R2, 1/R3,0,0,0; 0, (-(1/R2)-Kb), 1/R2, Kb, 0,0,0; 0, 1/R3,0,(-(1/R5)-(1/R3)-(1/R4)), (1/R5 +1/Zc), (1/R7),(-(1/R7)-(1/Zc)); 0,Kb,0,(-(1/R5)-Kb),(1/R5+1/Zc),0,-(1/Zc);0,0,0,0,0,(-(1/R6)-(1/R7)),1/R7;0,0,0,1,0,Kd/R6,-1];
Q = [1;0;0;0;0;0;0];
L = inv(G);
P = L*Q;
sprintf("%.6f\n",P)

Gain = abs(P(5));
Phase = angle(P(5));

t = 0:2e-6:20e-3;
v6f = Gain*sin(w*t+Phase);
f2 = figure();
plot(t*1000,v6f,"g");
xlabel("t[ms]");
ylabel("v6f[V]");
title("Forced Solution");
print(f2,"forced_solution.eps","-depsc");



fp5 = fopen('ngspice_4.txt',"w");
fprintf(fp5, "R1 1 2 %f\n",values(1));
fprintf(fp5, "R2 3 2 %f\n",values(2));
fprintf(fp5, "R3 2 5 %f\n",values(3));
fprintf(fp5, "R4 GND 5 %f\n",values(4));
fprintf(fp5, "R5 5 6 %f\n",values(5));
fprintf(fp5, "R6 GND 4 %f\n",values(6));
fprintf(fp5, "R7 7 8 %f\n",values(7));
fprintf(fp5, "C 6 8 %f\n",values(11));
fprintf(fp5, "Vs 1 GND sin(0 1 1k) \n");
fprintf(fp5, "V3 4 7 DC 0 \n");
fprintf(fp5, "Hc 5 8 V3 %f\n",values(10));
fprintf(fp5, "Gb 6 3 (2,5) %f\n",values(9));
fclose(fp5);

%total solution

t = -5e-3:2e-6:20e-3;
moment1 = t<0; moment2 = t>=0;
v6(moment1) = D(5);
v6(moment2) = X(4)*exp((-t(moment2))/Tau) + Gain*sin(w*(t(moment2))+Phase);
vs(moment1) = values(8);
vs(moment2) = sin(w*(t(moment2)));
f5 = figure();

plot(t*1000,v6,'k');
hold on
plot(t*1000,vs,'m');
xlabel("t[ms]");
 ylabel("vs_v6[V]");
 title("Total Solution");
legend('v6','vs','Location','Northeast');

print(f5,"Total_Solution.eps","-depsc");
	
	
	

