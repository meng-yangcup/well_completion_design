%参数设置井径，套管内径，套管外径,井深
D_h=444.5*10^-3;
d=313.6*10^-3;
D=339.7*10^-3;
h=260;
%水泥用量计算K1:井径扩大系数；rubber_stopper；pocket；
%
K1=1.1;
RS= 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huo水泥
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roucm = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;
roucm=1.82;rouw=1; m=0.515;Q=1.2;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;


%水利参数测算
%排量计算
fai_600=56;fai_300=40;
n=3.32*log10(fai_600-fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);


%计算分为3类
%雷诺数过渡流临界雷诺数；m/s;L/s
Re_c = 2100;
v_ci=(Re_c*800^n/(8000*roucm*d^n))^(1/(2-n));
v_ca=(0.2625*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
Q_ci=pi/40*d^2*v_ci;
Q_ca=pi/40*(D_h^2-D^2)*v_ca;
%塞流临界流速临界流速
Re_cp = 100;
v_ci_piston = (Re_cp*800^n*K/(8000*roucm*d^n))^(1/(2-n));
v_ca_piston = (0.0125*800^n*K/(roucm/(D_h^2-D^2)))^(1/(2-n));
Q_ca_piston = pi/40*(D_h^2-D^2)*v_ca_piston;
%套管鞋地层破裂压力rouc:钻井液密度h_cm水泥返高;h_套管井深，pf=MPa;
h_cm=0;
pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);
%计算施工最高泵压
p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
p_ha=10^-3*9.81*(roudf*(h-h_cm)+roucm*h_cm);
%净液柱压差
delta_p = p_ha - p_hi;

x = input('宾汉1Enter a number: ');

if x
    fi=0.03164/Re_c;
else
    a = (log10(n)+2.5)/50;
    b = (1.4-log10(n))/7;
    fi = a/Re^b;
end

%管内阻力
p_fi = 2*h*roudf*v_ci^2*fi/d;
%管外阻力
p_fa = 2*h*roudf*v_ca^2*fi/(D_h-D);

%summarise
p_max = p_ha +p_hi +p_fi +p_fa;


%%%注水泥时间
y = input('泵车123456Enter a number: ');

switch y
    case 1
        Q_pump = 19.9;
        p_rated = 33.1；
    case 2
        Q_pump = 23.4;
        p_rated = 28.2;
    case 3
        Q_pump = 27.1;
        p_rated = 24.3;
    case 4
        Q_pump = 31.1;
        p_rated = 21.2;
    case 5
        Q_pump = 35.4;
        p_rated = 18.6;
    case 6 
        Q_pump = 40;
        p_rated = 16.5;       
end
z = input('泵车数量Enter a number: ');
V_m = (pi/4)*d^2*h;
T_1 = V/(Q_pump*z);
T_2 = 3;
T_3 = V_m/Q_pump;
T_4 = 3;
T = T_1+ T_2+ T_3+ T_4;
fprintf('一开水泥返高:%f',h);
fprintf('一开水泥体积：%f',V);
fprintf('一开水泥浆密度:%f',roucm);
fprintf('一开干水泥重量kg：%f',W_c*1000);
fprintf('一开清水用量m^3:%f',Vw);
fprintf('一开石英砂用量kg：%f',W_c*1000/0.59*0.26);
fprintf('一开微珠用量kg：%f',W_c*1000/0.59*0.15);
fprintf('水灰比：%f',m);
%fprintf('水利参数设计结果：/n');
%fprintf('套缸直径mm：%f',);
%fprintf('额定冲次/min：%f',);
%fprintf('额定排量L/s：%f',);




%%%%%
%二开%
%二开%
%%%%%

%参数设置井径，套管内径，套管外径,井深
D_h=311.1*10^-3;
d=224.4*10^-3;
D=244.5*10^-3;
xx = input('班号Enter a number: ');
yy = input('学号后两位1234Enter a number: ');
h = 1500+(xx-1)*yy*3;

%水泥用量计算K1:井径扩大系数；rubber_stopper；pocket；
%
K1=1.1;
RS= 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huo水泥
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roucm = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;
roucm=1.82;rouw=1; m=0.515;Q=1.2;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;


%水利参数测算
%排量计算
fai_600=56;fai_300=40;
n=3.32*log10(fai_600-fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);


%计算分为3类
%雷诺数过渡流临界雷诺数；m/s;L/s
Re_c = 2100;
v_ci=(Re_c*800^n/(8000*roucm*d^n))^(1/(2-n));
v_ca=(0.2625*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
Q_ci=pi/40*d^2*v_ci;
Q_ca=pi/40*(D_h^2-D^2)*v_ca;
%塞流临界流速临界流速
Re_cp = 100;
v_ci_piston = (Re_cp*800^n*K/(8000*roucm*d^n))^(1/(2-n));
v_ca_piston = (0.0125*800^n*K/(roucm/(D_h^2-D^2)))^(1/(2-n));
Q_ca_piston = pi/40*(D_h^2-D^2)*v_ca_piston;
%套管鞋地层破裂压力rouc:钻井液密度h_cm水泥返高;h_套管井深，pf=MPa;
h_cm=0;
pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);
%计算施工最高泵压
p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
p_ha=10^-3*9.81*(roudf*(h-h_cm)+roucm*h_cm);
%净液柱压差
delta_p = p_ha - p_hi;

x = input('宾汉1Enter a number: ');

if x
    fi=0.03164/Re_c;
else
    a = (log10(n)+2.5)/50;
    b = (1.4-log10(n))/7;
    fi = a/Re^b;
end

%管内阻力
p_fi = 2*h*roudf*v_ci^2*fi/d;
%管外阻力
p_fa = 2*h*roudf*v_ca^2*fi/(D_h-D);

%summarise
p_max = p_ha +p_hi +p_fi +p_fa;


%%%注水泥时间
y = input('泵车1234Enter a number: ');

switch y
    case 1
        Q_pump = 19.9;
        p_rated = 33.1;
    case 2
        Q_pump = 23.4;
        p_rated = 28.2;
    case 3
        Q_pump = 27.1;
        p_rated = 24.3;
    case 4
        Q_pump = 31.1;
        p_rated = 21.2;
    case 5
        Q_pump = 35.4;
        p_rated = 18.6;
    case 6 
        Q_pump = 40;
        p_rated = 16.5;       
end
z = input('泵车数量Enter a number: ');
V_m = (pi/4)*d^2*h;
T_1 = V/(Q_pump*z);
T_2 = 3;
T_3 = V_m/Q_pump;
T_4 = 3;
T = T_1+ T_2+ T_3+ T_4;
fprintf('二开水泥返高:%f',h);
fprintf('二开水泥体积：%f',V);
fprintf('二开水泥浆密度:%f',roucm);
fprintf('二开干水泥重量kg：%f',W_c*1000);
fprintf('二开清水用量m^3:%f',Vw);
fprintf('二开石英砂用量kg：%f',W_c*1000/0.59*0.26);
fprintf('二开微珠用量kg：%f',W_c*1000/0.59*0.15);
fprintf('水灰比：%f',m);



%%%%%
%三开%
%三开%
%%%%%

%参数设置井径，套管内径，套管外径,井深
D_h=215.9*10^-3;
d=121.36*10^-3;
D=139.7*10^-3;
h = 2021;

%水泥用量计算K1:井径扩大系数；rubber_stopper；pocket；
%
K1=1.1;
RS= 4;
pocket=4;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huo水泥
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roucm = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;
roucm=1.82;rouw=1; m=0.515;Q=1.2;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;


%水利参数测算
%排量计算
fai_600=56;fai_300=40;
n=3.32*log10(fai_600-fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);


%计算分为3类
%雷诺数过渡流临界雷诺数；m/s;L/s
Re_c = 2100;
v_ci=(Re_c*800^n/(8000*roucm*d^n))^(1/(2-n));
v_ca=(0.2625*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
Q_ci=pi/40*d^2*v_ci;
Q_ca=pi/40*(D_h^2-D^2)*v_ca;
%塞流临界流速临界流速
Re_cp = 100;
v_ci_piston = (Re_cp*800^n*K/(8000*roucm*d^n))^(1/(2-n));
v_ca_piston = (0.0125*800^n*K/(roucm/(D_h^2-D^2)))^(1/(2-n));
Q_ca_piston = pi/40*(D_h^2-D^2)*v_ca_piston;
%套管鞋地层破裂压力rouc:钻井液密度h_cm水泥返高;h_套管井深，pf=MPa;
h_cm=0;
pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);
%计算施工最高泵压
p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
p_ha=10^-3*9.81*(roudf*(h-h_cm)+roucm*h_cm);
%净液柱压差
delta_p = p_ha - p_hi;

x = input('宾汉1Enter a number: ');

if x
    fi=0.03164/Re_c;
else
    a = (log10(n)+2.5)/50;
    b = (1.4-log10(n))/7;
    fi = a/Re^b;
end

%管内阻力
p_fi = 2*h*roudf*v_ci^2*fi/d;
%管外阻力
p_fa = 2*h*roudf*v_ca^2*fi/(D_h-D);

%summarise
p_max = p_ha +p_hi +p_fi +p_fa;


%%%注水泥时间
y = input('泵车1234Enter a number: ');

switch y
    case 1
        Q_pump = 19.9;
        p_rated = 33.1;
    case 2
        Q_pump = 23.4;
        p_rated = 28.2;
    case 3
        Q_pump = 27.1;
        p_rated = 24.3;
    case 4
        Q_pump = 31.1;
        p_rated = 21.2;
    case 5
        Q_pump = 35.4;
        p_rated = 18.6;
    case 6 
        Q_pump = 40;
        p_rated = 16.5;       
end
z = input('泵车数量Enter a number: ');
V_m = (pi/4)*d^2*h;
T_1 = V/(Q_pump*z);
T_2 = 3;
T_3 = V_m/Q_pump;
T_4 = 3;
T = T_1+ T_2+ T_3+ T_4;
fprintf('三开水泥返高:%f',h);
fprintf('三开水泥体积：%f',V);
fprintf('三开水泥浆密度:%f',roucm);
fprintf('三开干水泥重量kg：%f',W_c*1000);
fprintf('三开清水用量m^3:%f',Vw);
fprintf('三开石英砂用量kg：%f',W_c*1000/0.59*0.26);
fprintf('三开微珠用量kg：%f',W_c*1000/0.59*0.15);
fprintf('水灰比：%f',m);





