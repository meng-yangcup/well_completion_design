%参数设置井径，套管内径，套管外径,井深
D_h=444.5*10^-3;
d=313.6*10^-3;
D=339.7*10^-3;
h=260;
%水泥用量计算K1:井径扩大系数；rubber_stopper；pocket；
%
K1=1.1;
RS= 0;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%干水泥用量
K2=1.05;
rouc=3.16; rouw=1; m=0.515;Q=1.2;
q=rouc*rouw/(rouw+m*rouc)
W_c = K2* V *q
%水利参数测算