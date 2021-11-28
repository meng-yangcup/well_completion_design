%ahead
%参数设置井径，套管内径，套管外径,井深单位m
D_h=444.5*10^-3;
d=313.6*10^-3;
D=339.7*10^-3;
h=260;
%水泥用量计算K1:井径扩大系数； RS=rubber_stopper；pocket；
%
K1=1.1;
RS = 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huo水泥
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roudry = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;
roucm=1.82;rouw=1; m=0.515;
q=roudry*rouw/(rouw+m*roudry);
W_c = K2* V *q;   %%%t
Vw= m*W_c/rouw;   %%%m^3

%水利参数测算
%排量计算
fai_600=56;fai_300=40;
n=3.32*log10(fai_600/fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);
miudf_p = 4.6;


%%%选泵总类数量
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
    p_pump = (sum(c1(:,4)));%MPa%为简化计算认为泵压和排量可线性叠加
    Q_pump = (sum(c1(:,3)))/1000*60;%m^3/min假设泵效为100% %
    V_m = (pi/4)*d^2*h;% volume of mud

    T_1 = V/(Q_pump);
    T_2 = 3;
    T_3 = V_m/Q_pump;
    T_4 = 3;
    T = T_1+ T_2+ T_3+ T_4;

%计算分为2类紊流、断塞流
%%计算雷诺数：
%管内流动
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*(d*10)^n *vi^(2-n)/(800^n*K);

%环空流
    va = Q_pump*60/(pi/4*(D_h^2-D^2));
    Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
    %雷诺数过渡流临界雷诺数；m/s;L/s
        fprintf('水力参数合理，请选择流型');
        fprintf('环空流为紊流');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
%    v_ci=(Re_i*800^n/(8000*roudf*d^n))^(1/(2-n));%管柱内流速
%    v_ca=(Re_a*800^n*K/(8000*roucm/(D_h-D)^n))^(1/(2-n));%环空流速
%    Q_ci=pi/40*d^2*v_ci;                         %管柱排量
%    Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a < 100
        fprintf('环空流为塞流');
        fprintf('水力参数合理，请选择流型');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
    %v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
    %v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
    %Q_ci=pi/40*d^2*v_ci;
    %Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('水利参数不合理，请重新选泵');
        continue;
    end

%套管鞋地层破裂压力rouc:钻井液密度h_cm水泥返高;h_套管井深，pf=MPa;
    h_cm=0;
%pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:破裂密度
    rouf=100.0;
    pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
%计算施工最高泵压
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%净液柱压差
    delta_p = p_ha - p_hi;
    if Re_a > 2100
        x = input('宾汉1else幂律Enter a number: ');
        if x
            fi=0.03164/Re_a;
        else
            a = (log10(n)+2.5)/50;
            b = (1.4-log10(n))/7;
            fi = a/Re_a^b;
        end

%管内阻力
        p_fi = 2*h*roudf*1000*v_ci^2*fi/d;
%管外阻力
        p_fa = 2*h*roucm*1000*v_ca^2*fi/(D_h-D);
    else
        fi = 16/Re_i;
        fa = 24/Re_a;
    %%宾汉和幂律在环空流动计算很复杂；这么简化可能不对
%管内阻力
        p_fi = 2*h*roudf*1000*v_ci^2*fi/d;
%管外阻力
        p_fa = 2*h*roucm*1000*v_ca^2*fi/(D_h-D);
    end
%summarise
    p_max = delta_p +p_fi +p_fa;
    if p_pump - p_max < pf
        fprintf('泵车压力选择合理');
        break;
    else
        fprintf('压力参数不合理,请重新选泵');
    end
    
end

fprintf('一开水泥返高:%f\n',h);
fprintf('一开水泥体积：%f\n',V);
fprintf('一开水泥浆密度:%f\n',roucm);
fprintf('一开干水泥重量kg：%f\n',W_c*1000);
fprintf('一开清水用量m^3:%f\n',Vw);
fprintf('一开石英砂用量kg：%f\n',W_c*1000/0.59*0.26);
fprintf('一开微珠用量kg：%f\n',W_c*1000/0.59*0.15);
fprintf('水灰比：%f\n',m);
%fprintf('水利参数设计结果：/n');
for i = 1:3
    fprintf('套缸直径mm：%f',c1(i:1));
    fprintf('额定冲次/min：%f',c1(i,2));
    fprintf('额定排量L/s：%f',c1(i,3));
    fprintf('额定泵压MPa：%f',c1(i,4));
    fprintf('\n');
end
fprintf('注水泥施工时间：%f\n',T);
fprintf('替钻井液体积：%f\n',V_m);
fprintf('排量：%f\n',Q_pump);
fprintf('管内反速：%f\n',vi);
fprintf('环空反速：%f\n',va);

%%%%%
%二开%
%二开%
%%%%%


%参数设置井径，套管内径，套管外径,井深，重新赋值
D_h=311.1*10^-3;
d=224.4*10^-3;
D=244.5*10^-3;
xx = input('班号Enter a number: ');
yy = input('学号后两位1234Enter a number: ');
h = 1500+(xx-1)*50+yy*3;

%水泥用量计算K1:井径扩大系数；rubber_stopper； pocket；
%
K1=1.1;
RS= 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huo水泥
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roudry = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;
roucm=1.82;rouw=1; m=0.515;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;


%水利参数测算
%排量计算
%fai_600=56;fai_300=40;
%n=3.32*log10(fai_600-fai_300);
%K=0.511*fai_300/(511^n);
%miu_p = fai_600-fai_300;
%tao_0 = 0.511*(fai_300-miu_p);

%%%选泵总类数量
pf=0;p_pump=1;p_max=0;
while pf<p_pump-p_max
    c1=zeros(3,4);
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
    p_pump = (sum(c1(:,4)));%MPa
    Q_pump = (sum(c1(:,3)))/1000*60;%m^3/min
    V_m = (pi/4)*d^2*h;

    T_1 = V/(Q_pump);
    T_2 = 3;
    T_3 = V_m/Q_pump;
    T_4 = 3;
    T = T_1+ T_2+ T_3+ T_4;

%计算分为2类紊流、断塞流
%%计算雷诺数：
%管内流动
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*d^n *vi^(2-n)/(800^n*K);
%环空流
    va = Q_pump/60/(pi/4*(D_h^2-D^2));
    Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
    %雷诺数过渡流临界雷诺数；m/s;L/s
        fprintf('环空流为紊流');
        fprintf('水利参数合理，继续选择流型');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
        Q_ci=pi/40*d^2*v_ci;
        Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a < 100
        fprintf('环空流为塞流');
        fprintf('水利参数合理，继续选择流型');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
        Q_ci=pi/40*d^2*v_ci;
        Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('当前环空雷诺数%f ',Re_a);
        fprintf('水利参数不合理，请重新选泵');
        continue;
    end

%套管鞋地层破裂压力rouc:钻井液密度h_cm水泥返高;h_套管井深，pf=MPa;
    h_cm=0;
%pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:破裂密度
    rouf=10.0;
    pf=rouf* 9.81* (h-h_cm)/1000;%MPa
%计算施工最高泵压
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%净液柱压差
    delta_p = p_ha - p_hi;
    if Re_a > 2100

        x = input('宾汉1else幂律Enter a number: ');
        if x
            fi=0.03164/Re_a^0.25;
        else
            a = (log10(n)+2.5)/50;
            b = (1.4-log10(n))/7;
            fi = a/Re_a^b;
        end

    else
        fi = 16/Re_i;
        fa = 24/Re_a;
    %%宾汉和幂律在环空流动计算很复杂；这么简化可能不对
    end
%管内阻力
    p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;
%管外阻力
    p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;
%summarise
    p_max = delta_p +p_fi +p_fa;
    if p_pump-p_max < pf
        fprintf('压力参数合理');
        
    else
        fprintf('压力参数不合理，重选泵组');
    end
end
fprintf('二开水泥返高:%f\n',h);
fprintf('二开水泥体积：%f\n',V);
fprintf('二开水泥浆密度:%f\n',roucm);
fprintf('二开干水泥重量kg：%f\n',W_c*1000);
fprintf('二开清水用量m^3:%f\n',Vw);
fprintf('二开石英砂用量kg：%f\n',W_c*1000/0.59*0.26);
fprintf('二开微珠用量kg：%f\n',W_c*1000/0.59*0.15);
fprintf('水灰比：%f\n',m);

%fprintf('水利参数设计结果：/n');
for i = 1:3
    fprintf('套缸直径mm：%f',c1(i:1));
    fprintf('额定冲次/min：%f',c1(i,2));
    fprintf('额定排量L/s：%f',c1(i,3));
    fprintf('额定泵压MPa：%f',c1(i,4));
    fprintf('\n');
end
fprintf('注水泥施工时间：%f\n',T);
fprintf('替钻井液体积：%f\n',V_m);
fprintf('排量：%f\n',Q_pump);
fprintf('管内反速：%f\n',vi);
fprintf('环空反速：%f\n',va);

%%%%%
%三开%
%三开%
%%%%%

%参数设置井径，套管内径，套管外径,井深
D_h=215.9*10^-3;
d=121.36*10^-3;
D=139.7*10^-3;
h = 2021;

%水泥用量计算K1:井径扩大系数；rubber_stopper； pocket；
%
K1=1.1;
RS= 4;
pocket=4;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huo水泥
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;
roucm=1.82;rouw=1; m=0.515;Q=1.2;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;

%%%选泵总类数量
pf=0.1;p_pump=1;p_max=0;
while pf<p_pump-p_max
    c1=zeros(3,4);
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
    p_pump = (sum(c1(:,4)));%MPa
    Q_pump = (sum(c1(:,3)))/1000*60;%m^3/min
    V_m = (pi/4)*d^2*h;

    T_1 = V/(Q_pump);
    T_2 = 3;
    T_3 = V_m/Q_pump;
    T_4 = 3;
    T = T_1+ T_2+ T_3+ T_4;

%计算分为2类紊流、断塞流
%%计算雷诺数：
%管内流动
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*d^n *vi^(2-n)/(800^n*K);

%环空流
    va = Q_pump/60/(pi/4*(D_h^2-D^2));
    Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
    %雷诺数过渡流临界雷诺数；m/s;L/s
        fprintf('环空流为紊流');
        fprintf('继续选择流型');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
 %   Q_ci=pi/40*d^2*v_ci;
 %   Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a < 100
        fprintf('环空流为塞流');
        fprintf('继续选择流型');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
%    Q_ci=pi/40*d^2*v_ci;
%    Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('当前雷诺数为：%f',Re_a);
        fprintf('水利参数不合理，请重新选泵');
        continue;
    end

%套管鞋地层破裂压力rouc:钻井液密度h_cm水泥返高;h_套管井深，pf=MPa;
    h_cm=0;
    %pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:破裂密度
    rouf=2.0;
    pf=rouf* 9.81* (h-h_cm);
%计算施工最高泵压
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%净液柱压差
    delta_p = p_ha - p_hi;
    if Re_a > 2100

        x = input('宾汉1else幂律Enter a number: ');
        if x
            fi=0.03164/Re_a^0.25;
        else
           a = (log10(n)+2.5)/50;
           b = (1.4-log10(n))/7;
           fi = a/Re_a^b;
        end
    else
        fi = 16/Re_i;
        fa = 24/Re_a;
%%宾汉和幂律在环空流动计算很复杂；这么简化可能不对
%管内阻力
    p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;
%管外阻力
    p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;
    end
%summarise
    p_max = delta_p +p_fi +p_fa;

    if p_max - pf<pf
        fprintf('压力参数合理');
        break;
    else
        fprintf('压力参数不合理，重选泵组');
    end
end
fprintf('三开水泥返高:%f\n',h);
fprintf('三开水泥体积：%f\n',V);
fprintf('三开水泥浆密度:%f\n',roucm);
fprintf('三开干水泥重量kg：%f\n',W_c*1000);
fprintf('三开清水用量m^3:%f\n',Vw);
fprintf('三开石英砂用量kg：%f\n',W_c*1000/0.59*0.26);
fprintf('三开微珠用量kg：%f\n',W_c*1000/0.59*0.15);
fprintf('水灰比：%f\n',m);
%fprintf('水利参数设计结果：/n');
for i = 1:3
    fprintf('套缸直径mm：%f',c1(i:1));
    fprintf('额定冲次/min：%f',c1(i,2));
    fprintf('额定排量L/s：%f',c1(i,3));
    fprintf('额定泵压MPa：%f',c1(i,4));
    fprintf('\n');
end
fprintf('注水泥施工时间：%f\n',T);
fprintf('替钻井液体积：%f\n',V_m);
fprintf('排量：%f\n',Q_pump);
fprintf('管内反速：%f\n',vi);
fprintf('环空反速：%f\n',va);





