fprintf('一开/n');
D_h=444.5*10^-3;d=313.6*10^-3;D=339.7*10^-3;
h=260;K1=1.1;RS = 4;pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
V_m = (pi/4)*d^2*h;
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
roudry = rou_drycement;
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;roucm=1.82;rouw=1; m=0.515;
q=roudry*rouw/(rouw+m*roudry);
W_c = K2* V *q;   Vw= m*W_c/rouw;   %%%m^3

fai_600=56;fai_300=40;
n=3.32*log10(fai_600/fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);
miudf_p = 4.6;

c=[120,150,19.9,33.1,0;
   130,150,23.4,28.2,0;
   140,150,27.1,24.3,0;
   150,150,31.1,21.2,0;
   160,150,35.4,18.6,0;
   170,150,40.0,16.5,0;
   0 0 0 0 0];
an = [];
count=1;
opt = ones(3,5);
for i =1 :6
    for i1 = 1:7
        for i2 = 1:7
            c1=zeros(3,5);
            an(:,:,count)=[c(i,:);c(i1,:);c(i2,:)];
            p_pump = (sum(an(:,4,count)));
            Q_pump = (sum(an(:,3,count)))/1000*60;
            T_1 = V/(Q_pump);
            T_2 = 3;
            T_3 = V_m/Q_pump;
            T_4 = 3;
            T = T_1+ T_2+ T_3+ T_4;
            vi = Q_pump/60/(pi/4*d^2); %m/s
            Re_i = 8000*roudf*(d*10)^n *vi^(2-n)/(800^n*K);
            va = Q_pump/60/(pi/4*(D_h^2-D^2));
            Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
            h_cm=0;
            rouf=100.0;
            pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
            p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
            p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
            delta_p = p_ha - p_hi;
            fi=0.03164/Re_a;
            p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;    
            p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;    
            p_max = delta_p +p_fi +p_fa;

            x1=(pf-p_pump+p_max)/pf;
            an(:,5,count)=x1;
            c1(:,:,count)=an(:,:,count);

            if opt(:,5)>c1(:,5)
                opt=c1;
                fprintf('水泥返高:%f\n',h);
                fprintf('水泥体积：%f\n',V);
                fprintf('水泥浆密度:%f\n',roucm);
                fprintf('干水泥重量kg：%f\n',W_c*1000);
                fprintf('清水用量m^3:%f\n',Vw);
                fprintf('石英砂用量kg：%f\n',W_c*1000/0.59*0.26);
                fprintf('微珠用量kg：%f\n',W_c*1000/0.59*0.15);
                fprintf('水灰比：%f\n',m);
                fprintf('套缸直径mm:');
                fprintf('额定冲次/min：');
                fprintf('额定排量L/s：');
                fprintf('额定泵压MPa：');
                fprintf('\n');
                disp(opt);

                fprintf('注水泥施工时间：%f\n',T);
                fprintf('替钻井液体积：%f\n',V_m);
                fprintf('排量：%f\n',Q_pump);
                fprintf('管内反速：%f\n',vi);
                fprintf('环空反速：%f\n',va);
            end
            count=count+1;
            
        end
    end
end
fprintf('二开/n/n/n');
D_h=311.1*10^-3;d=224.4*10^-3;D=244.5*10^-3;
xx = input('班号Enter a number: ');
yy = input('学号后两位1234Enter a number: ');
h = 1500+(xx-1)*50+yy*3;

%水泥用量计算K1:井径扩大系数；rubber_stopper； pocket；
%
K1=1.1;RS= 4;pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huo水泥
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);roudry = rou_drycement;
rou_drillingfluid = 1.18;roudf=rou_drillingfluid;
%
K2=1.05;roucm=1.82;rouw=1; m=0.515;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;

fai_600=56;fai_300=40;
n=3.32*log10(fai_600/fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);
miudf_p = 4.6;

c=[120,150,19.9,33.1,0;
   130,150,23.4,28.2,0;
   140,150,27.1,24.3,0;
   150,150,31.1,21.2,0;
   160,150,35.4,18.6,0;
   170,150,40.0,16.5,0;
   0 0 0 0 0];
an = [];count=1;opt = ones(3,5);
for i =1 :6
    for i1 = 1:7
        for i2 = 1:7
            c1=zeros(3,5);
            an(:,:,count)=[c(i,:);c(i1,:);c(i2,:)];
            p_pump = (sum(an(:,4,count)));
            Q_pump = (sum(an(:,3,count)))/1000*60;
            T_1 = V/(Q_pump);
            T_2 = 3;
            T_3 = V_m/Q_pump;
            T_4 = 3;
            T = T_1+ T_2+ T_3+ T_4;
            vi = Q_pump/60/(pi/4*d^2); %m/s
            Re_i = 8000*roudf*(d*10)^n *vi^(2-n)/(800^n*K);
            va = Q_pump/60/(pi/4*(D_h^2-D^2));
            Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
            h_cm=0;
            rouf=100.0;
            pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
            p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
            p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
            delta_p = p_ha - p_hi;
            fi=0.03164/Re_a;
            p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;    
            p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;    
            p_max = delta_p +p_fi +p_fa;

            x1=(pf-p_pump+p_max)/pf;
            an(:,5,count)=x1;
            c1(:,:,count)=an(:,:,count);

            if opt(:,5)>c1(:,5)
                opt=c1;
                fprintf('水泥返高:%f\n',h);
                fprintf('水泥体积：%f\n',V);
                fprintf('水泥浆密度:%f\n',roucm);
                fprintf('干水泥重量kg：%f\n',W_c*1000);
                fprintf('清水用量m^3:%f\n',Vw);
                fprintf('石英砂用量kg：%f\n',W_c*1000/0.59*0.26);
                fprintf('微珠用量kg：%f\n',W_c*1000/0.59*0.15);
                fprintf('水灰比：%f\n',m);
                fprintf('套缸直径mm:');
                fprintf('额定冲次/min：');
                fprintf('额定排量L/s：');
                fprintf('额定泵压MPa：');
                fprintf('\n');
                disp(opt);

                fprintf('注水泥施工时间：%f\n',T);
                fprintf('替钻井液体积：%f\n',V_m);
                fprintf('排量：%f\n',Q_pump);
                fprintf('管内反速：%f\n',vi);
                fprintf('环空反速：%f\n',va);
            end
            count=count+1;
            
        end
    end
end
fprintf('三开/n/n/n');

D_h=215.9*10^-3;d=121.36*10^-3;D=139.7*10^-3;h = 2021;
K1=1.1;RS= 4;pocket=4;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
rou_drillingfluid = 1.18;roudf=rou_drillingfluid;
%
K2=1.05;roucm=1.82;rouw=1; m=0.515;Q=1.2;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;

fai_600=56;fai_300=40;
n=3.32*log10(fai_600/fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);
miudf_p = 4.6;

c=[120,150,19.9,33.1,0;
   130,150,23.4,28.2,0;
   140,150,27.1,24.3,0;
   150,150,31.1,21.2,0;
   160,150,35.4,18.6,0;
   170,150,40.0,16.5,0;
   0 0 0 0 0];
an = [];count=1;opt = ones(3,5);
for i =1 :6
    for i1 = 1:7
        for i2 = 1:7
            c1=zeros(3,5);
            an(:,:,count)=[c(i,:);c(i1,:);c(i2,:)];
            p_pump = (sum(an(:,4,count)));
            Q_pump = (sum(an(:,3,count)))/1000*60;
            T_1 = V/(Q_pump);
            T_2 = 3;
            T_3 = V_m/Q_pump;
            T_4 = 3;
            T = T_1+ T_2+ T_3+ T_4;
            vi = Q_pump/60/(pi/4*d^2); %m/s
            Re_i = 8000*roudf*(d*10)^n *vi^(2-n)/(800^n*K);
            va = Q_pump/60/(pi/4*(D_h^2-D^2));
            Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
            h_cm=0;
            rouf=100.0;
            pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
            p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
            p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
            delta_p = p_ha - p_hi;
            fi=0.03164/Re_a;
            p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;    
            p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;    
            p_max = delta_p +p_fi +p_fa;

            x1=(pf-p_pump+p_max)/pf;
            an(:,5,count)=x1;
            c1(:,:,count)=an(:,:,count);

            if opt(:,5)>c1(:,5)
                opt=c1;
                fprintf('水泥返高:%f\n',h);
                fprintf('水泥体积：%f\n',V);
                fprintf('水泥浆密度:%f\n',roucm);
                fprintf('干水泥重量kg：%f\n',W_c*1000);
                fprintf('清水用量m^3:%f\n',Vw);
                fprintf('石英砂用量kg：%f\n',W_c*1000/0.59*0.26);
                fprintf('微珠用量kg：%f\n',W_c*1000/0.59*0.15);
                fprintf('水灰比：%f\n',m);
                fprintf('套缸直径mm:');
                fprintf('额定冲次/min：');
                fprintf('额定排量L/s：');
                fprintf('额定泵压MPa：');
                fprintf('\n');
                disp(opt);

                fprintf('注水泥施工时间：%f\n',T);
                fprintf('替钻井液体积：%f\n',V_m);
                fprintf('排量：%f\n',Q_pump);
                fprintf('管内反速：%f\n',vi);
                fprintf('环空反速：%f\n',va);
            end
            count=count+1;
            
        end
    end
end


