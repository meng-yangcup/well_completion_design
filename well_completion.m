%Set initial value;
D_h=444.5*10^-3; %D_h:Well diameter:unit m;
d=313.6*10^-3;   %d:Inner diameter of casing:unit m;
D=339.7*10^-3;   %D:Outer diameter of casing:unit m;
h=260;           %h:depth of well:unit m
%Calculate the amount of cement unit m^3;
%K1:Diameter expansion factor; RS=rubber_stopperunit m;pocketunit m;
K1=1.1; RS = 4; pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%Weighted average calculation of dry cement density;
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roudry = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%Ground loss factor;
K2=1.05;
roucm=1.82;rouw=1; 
%syms V_w V_c Volume_of_Water;Volume_of_Cement;
%[V_w,V_c]=solve('V_w + V_c = V','V_w + 2.65*V_c = V*1.82');
A=[1 1;1 2.65];B=[V;1.82*V];
X=A\B;
V_w=X(1);
V_c=X(2);
m=V_w/2.65/V_c;  %Water-cement ratio;
q=roudry*rouw/(rouw+m*roudry);
W_c = K2* V *q;   %%%t
Vw= m*W_c/rouw;   %%%m^3

% Water conservancy parameter calculation
% Displacement calculation
fai_600=56;fai_300=40;
n=3.32*log10(fai_600/fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);
miudf_p = 4.6;


%%%Select the number and the types of pumps
c=[120,150,19.9,33.1;
   130,150,23.4,28.2;
   140,150,27.1,24.3;
   150,150,31.1,21.2;
   160,150,35.4,18.6;
   170,150,40.0,16.5];
c1=zeros(3,4);
pf=0;p_pump=1;p_max=0;
while pf<p_pump-p_max
    for i=1:3
        k = input('泵车0123456Enter a number: ');
        if k==0
            continue
        else
            for j=i:3
                c1(j,:)=c(k,:);
            end
        end
    end
    %In order to simplify the calculation, 
    %it is considered that the pump pressure and 
    %displacement can be linearly superimposed;
    %Assuming that the pump efficiency is 100%;
    p_pump = (sum(c1(:,4)));%MPa
    Q_pump = (sum(c1(:,3)))/1000*60; %m^3/min 
    V_m = (pi/4)*d^2*h;              %volume_of_mud
    T_1 = V/(Q_pump);
    T_2 = 3;
    T_3 = V_m/Q_pump;
    T_4 = 3;
    T = T_1+ T_2+ T_3+ T_4;

%The calculation is divided into 2 types of
%turbulent flow and blocked flow;
%%Calculate the Reynolds number:
%Flow in the Tube
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*(100*d)^n *vi^(2-n)/(800^n*K);

%Flow in the Annulus 
    va = Q_pump/60/(pi/4*(D_h^2-D^2));%m/s
    Re_a = 8000*roucm*(100*(D_h-D))^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
%Critical Reynolds number of transition flow；m/s;
        fprintf('水力参数合理，请选择流型/n');
        fprintf('环空流为紊流');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
        if v_ca>2
            fprintf('环空反速过快，/n');
            continue;
        end
%    Another way of thinking:
%    Assuming that the Re_i=2100;Calculate the Critical velocity
%    Critical quantity of flow; choose the pumps that meet the requirments;
%    v_ci=(Re_i*800^n/(8000*roudf*d^n))^(1/(2-n));%
%    v_ca=(Re_a*800^n*K/(8000*roucm/(D_h-D)^n))^(1/(2-n));
%    Q_ci=pi/40*d^2*v_ci;
%    Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a <500%Re<500 /~=0
        fprintf('环空流为层流/塞流, ''/n');
        fprintf('水力参数合理，请选择流型');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
    %v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
    %v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
    %Q_ci=pi/40*d^2*v_ci;
    %Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('水利参数不合理，请重新选泵');
        fprintf('当前雷诺数%f\n',Re_a);
        continue;
    end

%rouc:drilling fluid density; h_cm:Cement anti-high;h:Casing depth
    h_cm=0;
%pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:density of fracture
    rouf=100.0;roul=1.0;
    pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
    pl=roul*9.8*(h-h_cm)*0.001;
%Calculate the highest pump pressure for construction
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%Net liquid column pressure difference
    delta_p = p_ha - p_hi;
    if delta_p>pl
        fprintf('存在漏失风险，');
        continue;
    else
        fprintf('无漏失风险，');
    end
    if Re_a > 2100
        x = input('宾汉1else幂律Enter a number: ');
        if x
            fi=0.03164/(Re_a^0.25);
        else
            a = (log10(n)+2.5)/50;
            b = (1.4-log10(n))/7;
            fi = a/Re_a^b;
        end
%Resistance in the tube;
        p_fi = 2*h*roudf*vi^2*fi/d*10^-3;
%Resistance in the annulus;
        p_fa = 2*h*roucm*va^2*fi/(D_h-D)*10^-3;
    else
        fi = 16/Re_i;
        fa = 24/Re_a;
%%The calculation of Bingham and power
%law flowing in the annulus is very complicated; 
%this simplification may be wrong
%Resistance in the tube;
        p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;
%Resistance in the annulus;
        p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;
    end
%summarise
    p_max = delta_p +p_fi +p_fa;
    p_pumpneed=pf + p_max;
    if p_pump >p_pumpneed%%%一开不考虑压裂地层
        fprintf('泵车压力足够/n');
        
    else
        fprintf('压力参数不合理,压裂地层，请重新选泵');
        continue;
    end
    answer1 = {'水泥返高',h 
    '水泥体积:',V;
    '水泥浆密度:',roucm;
    '干水泥重量kg:',W_c*1000;
    '干水泥重量kg:',W_c*1000;
    '清水用量m^3:',Vw;
    '石英砂用量kg:',W_c*1000/0.59*0.26;
    '微珠用量kg:',W_c*1000/0.59*0.15;
    '水灰比:',m;
    '套缸直径mm:',c1(:,1);
    '额定冲次/min:',c1(:,2);
    '额定排量L/s:',c1(:,3);
    '额定泵压MPa:',c1(:,4);
    '注水泥施工时间:',T;
    '替钻井液体积:',V_m;
    '排量:',Q_pump;
    '管内反速:',vi;
    '环空反速:',va;
    '所需泵压',p_pumpneed;
    '雷诺数',Re_a};
end

fprintf('一开水泥返高:%f\n',h);
fprintf('一开水泥体积：%f\n',V);
fprintf('一开水泥浆密度:%f\n',roucm);
fprintf('一开干水泥重量kg：%f\n',W_c*1000*0.59);
fprintf('一开清水用量m^3:%f\n',Vw);
fprintf('一开石英砂用量kg：%f\n',W_c*1000*0.26);
fprintf('一开微珠用量kg：%f\n',W_c*1000*0.15);
fprintf('水灰比：%f\n',m);
%fprintf('水利参数设计结果：/n');
for i = 1:3
    fprintf('套缸直径mm：%f',c1(i:1));
    fprintf('额定冲次/min：%f',c1(i,2));
    fprintf('额定排量L/s：%f',c1(i,3));
    fprintf('额定泵压MPa：%f',c1(i,4));
    fprintf('\n');
end
fprintf('当前雷诺数%f\n',Re_a);
fprintf('注水泥施工时间：%f\n',T);
fprintf('替钻井液体积：%f\n',V_m);
fprintf('排量：%f\n',Q_pump);
fprintf('管内反速：%f\n',vi);
fprintf('环空反速：%f\n',va);
fprintf('雷诺数%f',Re_a);

%%%%%
%二开%
%二开%
%%%%%
%Parameter setting well diameter, casing inner diameter,
%casing outer diameter, well depth, re-assignment
D_h=311.1*10^-3;
d=224.4*10^-3;
D=244.5*10^-3;
xx = input('班号Enter a number: ');
yy = input('学号后两位1234Enter a number: ');
h = 1500+(xx-1)*50+yy*3;

%Cement consumption calculation
%K1:Diameter expansion factor;rubber_stopper； pocket;
K1=1.1;
RS= 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%Calculate the amount of cement unit m^3;
%K1:Diameter expansion factor; RS=rubber_stopperunit m;pocketunit m;
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roudry = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%%Ground loss factor;
K2=1.05;
A=[1 1;1 2.65];B=[V;1.82*V];
X=A\B;
V_w=X(1);
V_c=X(2);
m=V_w/2.65/V_c;  %Water-cement ratio;
q=roudry*rouw/(rouw+m*roudry);
W_c = K2* V *q;   %%%t
Vw= m*W_c/rouw;


%%%Select the number and the types of pumps
c=[120,150,19.9,33.1;
   130,150,23.4,28.2;
   140,150,27.1,24.3;
   150,150,31.1,21.2;
   160,150,35.4,18.6;
   170,150,40.0,16.5];
c1=zeros(3,4);
pf=0;p_pump=1;p_max=0;
while pf<p_pump-p_max
    for i=1:3
        k = input('泵车0123456Enter a number: ');
        if k==0
            continue
        else
            for j=i:3
                c1(j,:)=c(k,:);
            end
        end
    end
    %In order to simplify the calculation, 
    %it is considered that the pump pressure and 
    %displacement can be linearly superimposed;
    %Assuming that the pump efficiency is 100%;
    p_pump = (sum(c1(:,4)));%MPa
    Q_pump = (sum(c1(:,3)))/1000*60; %m^3/min 
    V_m = (pi/4)*d^2*h;              %volume_of_mud
    T_1 = V/(Q_pump);
    T_2 = 3;
    T_3 = V_m/Q_pump;
    T_4 = 3;
    T = T_1+ T_2+ T_3+ T_4;

%The calculation is divided into 2 types of
%turbulent flow and blocked flow;
%%Calculate the Reynolds number:
%Flow in the Tube
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*(100*d)^n *vi^(2-n)/(800^n*K);

%Flow in the Annulus 
    va = Q_pump/60/(pi/4*(D_h^2-D^2));%m/s
    Re_a = 8000*roucm*(100*(D_h-D))^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
%Critical Reynolds number of transition flow；m/s;
        fprintf('水力参数合理，请选择流型/nRe_a%f',Re_a);
        fprintf('环空流为紊流');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
        if va>2
            fprintf('环空反速过快，/n');
            continue;
        end
%    Another way of thinking:
%    Assuming that the Re_i=2100;Calculate the Critical velocity
%    Critical quantity of flow; choose the pumps that meet the requirments;
%    v_ci=(Re_i*800^n/(8000*roudf*d^n))^(1/(2-n));%
%    v_ca=(Re_a*800^n*K/(8000*roucm/(D_h-D)^n))^(1/(2-n));
%    Q_ci=pi/40*d^2*v_ci;
%    Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a <500%Re<500 /~=0
        fprintf('环空流为层流/塞流, ''/n');
        fprintf('水力参数合理，请选择流型');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
    %v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
    %v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
    %Q_ci=pi/40*d^2*v_ci;
    %Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('水利参数不合理，请重新选泵');
        fprintf('当前雷诺数%f\n',Re_a);
        continue;
    end

%rouc:drilling fluid density; h_cm:Cement anti-high;h:Casing depth
    h_cm=0;
%pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:density of fracture
    rouf=2.0;roul=1.0;
    pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
    pl=roul*9.8*(h-h_cm)*0.001;
%Calculate the highest pump pressure for construction
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%Net liquid column pressure difference
    delta_p = p_ha - p_hi;
    if delta_p>pl
        fprintf('存在漏失风险，');
        continue;
    else
        fprintf('无漏失风险，');
    end
    if Re_a > 2100
        x = input('宾汉1else幂律Enter a number: ');
        if x
            fi=0.03164/(Re_a^0.25);
        else
            a = (log10(n)+2.5)/50;
            b = (1.4-log10(n))/7;
            fi = a/Re_a^b;
        end
%Resistance in the tube;
        p_fi = 2*h*roudf*vi^2*fi/d*10^-3;
%Resistance in the annulus;
        p_fa = 2*h*roucm*va^2*fi/(D_h-D)*10^-3;
    else
        fi = 16/Re_i;
        fa = 24/Re_a;
%%The calculation of Bingham and power
%law flowing in the annulus is very complicated; 
%this simplification may be wrong
%Resistance in the tube;
        p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;
%Resistance in the annulus;
        p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;
    end
%summarise
    p_max = delta_p +p_fi +p_fa;
    p_pumpneed=pf + p_max;
    if p_pump < p_pumpneed
        fprintf('泵车压力足够/n');
        
    else
        fprintf('压力参数不合理,压裂地层，请重新选泵');
        continue;
    end
    answer2 = {'水泥返高',h 
    '水泥体积:',V;
    '水泥浆密度:',roucm;
    '干水泥重量kg:',W_c*1000*0.59;
    '清水用量m^3:',Vw;
    '石英砂用量kg:',W_c*1000*0.26;
    '微珠用量kg:',W_c*1000*0.15;
    '水灰比:',m;
    '套缸直径mm:',c1(:,1);
    '额定冲次/min:',c1(:,2);
    '额定排量L/s:',c1(:,3);
    '额定泵压MPa:',c1(:,4);
    '注水泥施工时间:',T;
    '替钻井液体积:',V_m;
    '排量:',Q_pump;
    '管内反速:',vi;
    '环空反速:',va;
    '所需泵压',p_pumpneed;
    '雷诺数',Re_a};
end
fprintf('二开水泥返高:%f\n',h);
fprintf('二开水泥体积：%f\n',V);
fprintf('二开水泥浆密度:%f\n',roucm);
fprintf('二开干水泥重量kg：%f\n',W_c*1000*0.59);
fprintf('二开清水用量m^3:%f\n',Vw);
fprintf('二开石英砂用量kg：%f\n',W_c*1000*0.26);
fprintf('二开微珠用量kg：%f\n',W_c*1000*0.15);
fprintf('水灰比：%f\n',m);
fprintf('雷诺数%f',Re_a);
for i = 1:3
    fprintf('套缸直径mm：%f',c1(i:1));
    fprintf('额定冲次/min：%f',c1(i,2));
    fprintf('额定排量L/s：%f',c1(i,3));
    fprintf('额定泵压MPa：%f',c1(i,4));
    fprintf('\n');
end
fprintf('当前雷诺数%f\n',Re_a);
fprintf('注水泥施工时间：%f\n',T);
fprintf('替钻井液体积：%f\n',V_m);
fprintf('排量：%f\n',Q_pump);
fprintf('管内反速：%f\n',vi);
fprintf('环空反速：%f\n',va);

%%%%%
%三开%
%三开%
%%%%%
D_h=215.9*10^-3;
d=121.36*10^-3;
D=139.7*10^-3;
h = 2021;
K1=1.1;
RS= 4;
pocket=4;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
K2=1.05;
A=[1 1;1 2.65];B=[V;1.82*V];
X=A\B;
V_w=X(1);
V_c=X(2);
m=V_w/2.65/V_c; 
q=roudry*rouw/(rouw+m*roudry);
W_c = K2* V *q; 
Vw= m*W_c/rouw;


%%%Select the number and the types of pumps
c1=zeros(3,4);
pf=0;p_pump=1;p_max=0;
while pf<p_pump-p_max
c=[120,150,19.9,33.1;
   130,150,23.4,28.2;
   140,150,27.1,24.3;
   150,150,31.1,21.2;
   160,150,35.4,18.6;
   170,150,40.0,16.5];
    for i=1:3
        k = input('泵车0123456Enter a number: ');
        if k==0
            continue
        else
            for j=i:3
                c1(j,:)=c(k,:);
            end
        end
    end
    %In order to simplify the calculation, 
    %it is considered that the pump pressure and 
    %displacement can be linearly superimposed;
    %Assuming that the pump efficiency is 100%;
    p_pump = (sum(c1(:,4)));%MPa
    Q_pump = (sum(c1(:,3)))/1000*60; %m^3/min 
    V_m = (pi/4)*d^2*h;              %volume_of_mud
    T_1 = V/(Q_pump);
    T_2 = 3;
    T_3 = V_m/Q_pump;
    T_4 = 3;
    T = T_1+ T_2+ T_3+ T_4;

%The calculation is divided into 2 types of
%turbulent flow and blocked flow;
%%Calculate the Reynolds number:
%Flow in the Tube
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*(100*d)^n *vi^(2-n)/(800^n*K);

%Flow in the Annulus 
    va = Q_pump/60/(pi/4*(D_h^2-D^2));%m/s
    Re_a = 8000*roucm*(100*(D_h-D))^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
%Critical Reynolds number of transition flow；m/s;
        fprintf('水力参数合理，请选择流型/n');
        fprintf('环空流为紊流');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
        if va>2
            fprintf('环空反速过快，/n');
            continue;
        end
%    Another way of thinking:
%    Assuming that the Re_i=2100;Calculate the Critical velocity
%    Critical quantity of flow; choose the pumps that meet the requirments;
%    v_ci=(Re_i*800^n/(8000*roudf*d^n))^(1/(2-n));%
%    v_ca=(Re_a*800^n*K/(8000*roucm/(D_h-D)^n))^(1/(2-n));
%    Q_ci=pi/40*d^2*v_ci;
%    Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a <500%Re<500 /~=0
        fprintf('环空流为层流/塞流, ''/n');
        fprintf('水力参数合理，请选择流型');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
    %v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
    %v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
    %Q_ci=pi/40*d^2*v_ci;
    %Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('水利参数不合理，请重新选泵');
        fprintf('当前雷诺数%f\n',Re_a);
        continue;
    end

%rouc:drilling fluid density; h_cm:Cement anti-high;h:Casing depth
    h_cm=0;
%pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:density of fracture
    rouf=2.0;roul=1.0;
    pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
    pl=roul*9.8*(h-h_cm)*0.001;
%Calculate the highest pump pressure for construction
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%Net liquid column pressure difference
    delta_p = p_ha - p_hi;
    if delta_p>pl
        fprintf('存在漏失风险，');
        continue;
    else
        fprintf('无漏失风险，');
    end
    if Re_a > 2100
        x = input('宾汉1else幂律Enter a number: ');
        if x
            fi=0.03164/(Re_a^0.25);
        else
            a = (log10(n)+2.5)/50;
            b = (1.4-log10(n))/7;
            fi = a/Re_a^b;
        end
%Resistance in the tube;
        p_fi = 2*h*roudf*vi^2*fi/d*10^-3;
%Resistance in the annulus;
        p_fa = 2*h*roucm*va^2*fi/(D_h-D)*10^-3;
    else
        fi = 16/Re_i;
        fa = 24/Re_a;
%%The calculation of Bingham and power
%law flowing in the annulus is very complicated; 
%this simplification may be wrong
%Resistance in the tube;
        p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;
%Resistance in the annulus;
        p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;
    end
%summarise
    p_max = delta_p +p_fi +p_fa;
    p_pumpneed=pf + p_max;
    if p_pump < p_pumpneed
        fprintf('泵车压力足够/n');
        
    else
        fprintf('压力参数不合理,压裂地层，请重新选泵');
        continue;
    end
    answer3 = {'水泥返高',h 
    '水泥体积:',V;
    '水泥浆密度:',roucm;
    '干水泥重量kg:',W_c*1000;
    '干水泥重量kg:',W_c*1000*0.59;
    '清水用量m^3:',Vw;
    '石英砂用量kg:',W_c*1000*0.26;
    '微珠用量kg:',W_c*1000*0.15;
    '水灰比:',m;
    '套缸直径mm:',c1(:,1);
    '额定冲次/min:',c1(:,2);
    '额定排量L/s:',c1(:,3);
    '额定泵压MPa:',c1(:,4);
    '注水泥施工时间:',T;
    '替钻井液体积:',V_m;
    '排量:',Q_pump;
    '管内反速:',vi;
    '环空反速:',va;
    '所需泵压',p_pumpneed
    '雷诺数',Re_a};
end
fprintf('三开水泥返高:%f\n',h);
fprintf('三开水泥体积：%f\n',V);
fprintf('三开水泥浆密度:%f\n',roucm);
fprintf('三开干水泥重量kg：%f\n',W_c*1000*0.59);
fprintf('三开清水用量m^3:%f\n',Vw);
fprintf('三开石英砂用量kg：%f\n',W_c*1000*0.26);
fprintf('三开微珠用量kg：%f\n',W_c*1000*0.15);
fprintf('水灰比：%f\n',m);
%fprintf('水利参数设计结果：/n');
for i = 1:3
    fprintf('套缸直径mm：%f',c1(i:1));
    fprintf('额定冲次/min：%f',c1(i,2));
    fprintf('额定排量L/s：%f',c1(i,3));
    fprintf('额定泵压MPa：%f',c1(i,4));
    fprintf('\n');
end
fprintf('当前雷诺数%f\n',Re_a);
fprintf('注水泥施工时间：%f\n',T);
fprintf('替钻井液体积：%f\n',V_m);
fprintf('排量：%f\n',Q_pump);
fprintf('管内反速：%f\n',vi);
fprintf('环空反速：%f\n',va);





