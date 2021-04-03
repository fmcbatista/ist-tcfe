%Nós
R1 = 1.04001336091 *10^3
R2 = 2.04372276851 *10^3
R3 = 3.11359737601 *10^3
R4 = 4.17085404861 *10^3
R5 = 3.02859283303 *10^3
R6 = 2.070545767 *10^3
R7 = 1.01835949725 *10^3
Va = 5.20102702949 
Id = 1.00460501759 *10^-3
Kb = 7.19043597753 *10^-3
Kc = 8.06397385506*10^3


A=[0;0;-Va;Id;Id;0;0]


B=[(-1/R2)-(1/R3)-(1/R1),1/R2,0,1/R3,0,0,0;Kb+(1/R2),-1/R2,0,-Kb,0,0,0;0,0,1,0,0,0,0;1/R3,0,1/R4,(-1/R4)+(-1/R3)-(1/R5),1/R5,1/R7,-1/R7;Kb,0,0,(-1/R5)-(Kb),1/R5,0,0;0,0,1/R6,0,0,(-1/R6)-(1/R7),1/R7;0,0,Kc/R6,-1,0,-Kc/R6,1]

C=inv(B)

D=C*A

v1=D(1);
v2=D(2);
v3=D(3);
v4=D(4);
v5=D(5);
v6=D(6);
v7=D(7);

fp=fopen('/home/fmcb/ist-tcfe/t1/doc/ValoresNos.tex',"w");
fprintf(fp,"v1 & %f \\\\ \\hline\nv2 & %f \\\\ \\hline\nv3 & %f \\\\ \\hline\nv4 & %f \\\\ \\hline\nv5 & %f \\\\ \\hline\nv6 & %f \\\\ \\hline\nv7 & %f \\\\ \\hline",v1,v2,v3,v4,v5,v6,v7);
fclose(fp);

%Malhas

E=[R1+R2+R3,R3,R4;R4,0,R7+R6+R4-Kc;R3,R3-(1/Kb),0]

F=[Va;0;0]

G=inv(E)

H=G*F

i1=H(1);
i2=H(2);
i3=H(3);

fq=fopen('/home/fmcb/ist-tcfe/t1/doc/ValoresMalhas.tex',"w");
fprintf(fq,"Ia & %f \\\\ \\hline\nIb & %f \\\\ \\hline\nIc & %f \\\\ \\hline",i1,i2,i3);
fclose(fq);
