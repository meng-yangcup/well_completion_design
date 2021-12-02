%ahead
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

c=[120,150,19.9,33.1,0;
   130,150,23.4,28.2,0;
   140,150,27.1,24.3,0;
   150,150,31.1,21.2,0;
   160,150,35.4,18.6,0;
   170,150,40.0,16.5,0;
   0 0 0 0 0];
answer = zeros(20,5);
an = [];
count=1;
opt = ones(3,5);
for i0 =1 :6
    for i1 = 1:7
        for i2 = 1:7
            c1=zeros(3,5);
            an(:,:,count)=[c(i0,:);c(i1,:);c(i2,:)];
            p_pump = (sum(an(:,4,count)));
            Q_pump = (sum(an(:,3,count)))/1000*60;

            V_m = (pi/4)*d^2*h;              %volume_of_mud
            T_1 = V/(Q_pump);
            T_2 = 3;
            T_3 = V_m/Q_pump;
            T_4 = 3;
            T = T_1+ T_2+ T_3+ T_4;

            vi = Q_pump/60/(pi/4*d^2); %m/s
            Re_i = 8000*roudf*(100*d)^n *vi^(2-n)/(800^n*K);

            va = Q_pump/60/(pi/4*(D_h^2-D^2));%m/s
            Re_a = 8000*roucm*(100*(D_h-D))^n *va^(2-n)/(800^n*K);
            if Re_a > 2100
                v_ci = Q_pump/(pi/4*d^2);
                v_ca = Q_pump/(pi/4*(D_h^2-D^2));
            else %Re_a <500%Re<500 /~=0 elseif
                v_ci = Q_pump/(pi/4*d^2);
                v_ca = Q_pump/(pi/4*(D_h^2-D^2));
%            else
%                continue;
            end
            %rouc:drilling fluid density; h_cm:Cement anti-high;h:Casing depth
            h_cm=0;
        %pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:density of fracture
            rouf=10.0;
            pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
        %Calculate the highest pump pressure for construction
            p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
            p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
        %Net liquid column pressure difference
            delta_p = p_ha - p_hi;
            if Re_a > 2100
                x = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%1幂律2宾汉
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
                p_fi = 2*h*roudf*1000*v_ci^2*fi/d;
        %Resistance in the annulus;
                p_fa = 2*h*roucm*1000*v_ca^2*fi/(D_h-D);
            end
        %summarise
            p_max = delta_p +p_fi +p_fa;
            if p_pump - p_max < pf
                x1=(pf-p_pump+p_max)/pf;
                an(:,5,count)=x1;
                c1(:,:,count)=an(:,:,count);
                if opt(:,5)>c1(:,5)
                    opt=c1;
                end
%            else
%            continue;
            end

            x1=(pf-p_pump+p_max)/pf;
            an(:,5,count)=x1;
            c1(:,:,count)=an(:,:,count);
            if opt(:,5)>c1(:,5)
                opt=c1;
            end
            count=count+1;
            answer = {'水泥返高',h;
                '水泥体积:',V;
                '水泥浆密度:',roucm;
                '干水泥重量kg:',W_c*1000;
                '干水泥重量kg:',W_c*1000;
                '清水用量m^3:',Vw;
                '石英砂用量kg:',W_c*1000/0.59*0.26;
                '微珠用量kg:',W_c*1000/0.59*0.15;
                '水灰比:',m;
                '套缸直径mm:',an(:,1,count-1);
                '额定冲次/min:',an(:,2,count-1);
                '额定排量L/s:',an(:,3,count-1);
                '额定泵压MPa:',an(:,4,count-1);
                '注水泥施工时间:',T;
                '替钻井液体积:',V_m;
                '排量:',Q_pump;
                '管内反速:',vi;
                '环空反速:',va};
        end
    end
end
