%�������þ������׹��ھ����׹��⾶,����
D_h=444.5*10^-3;
d=313.6*10^-3;
D=339.7*10^-3;
h=260;
%ˮ����������K1:��������ϵ����rubber_stopper��pocket��
%
K1=1.1;
RS= 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huoˮ��
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


%ˮ����������
%��������
fai_600=56;fai_300=40;
n=3.32*log10(fai_600-fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);


%�����Ϊ3��
%��ŵ���������ٽ���ŵ����m/s;L/s
Re_c = 2100;
v_ci=(Re_c*800^n/(8000*roucm*d^n))^(1/(2-n));
v_ca=(0.2625*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
Q_ci=pi/40*d^2*v_ci;
Q_ca=pi/40*(D_h^2-D^2)*v_ca;
%�����ٽ������ٽ�����
Re_cp = 100;
v_ci_piston = (Re_cp*800^n*K/(8000*roucm*d^n))^(1/(2-n));
v_ca_piston = (0.0125*800^n*K/(roucm/(D_h^2-D^2)))^(1/(2-n));
Q_ca_piston = pi/40*(D_h^2-D^2)*v_ca_piston;
%�׹�Ь�ز�����ѹ��rouc:�꾮Һ�ܶ�h_cmˮ�෵��;h_�׹ܾ��pf=MPa;
h_cm=0;
pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);
%����ʩ����߱�ѹ
p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
p_ha=10^-3*9.81*(roudf*(h-h_cm)+roucm*h_cm);
%��Һ��ѹ��
delta_p = p_ha - p_hi;

x = input('����1Enter a number: ');

if x
    fi=0.03164/Re_c;
else
    a = (log10(n)+2.5)/50;
    b = (1.4-log10(n))/7;
    fi = a/Re^b;
end

%��������
p_fi = 2*h*roudf*v_ci^2*fi/d;
%��������
p_fa = 2*h*roudf*v_ca^2*fi/(D_h-D);

%summarise
p_max = p_ha +p_hi +p_fi +p_fa;


%%%עˮ��ʱ��
y = input('�ó�123456Enter a number: ');

switch y
    case 1
        Q_pump = 19.9;
        p_rated = 33.1��
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
z = input('�ó�����Enter a number: ');
V_m = (pi/4)*d^2*h;
T_1 = V/(Q_pump*z);
T_2 = 3;
T_3 = V_m/Q_pump;
T_4 = 3;
T = T_1+ T_2+ T_3+ T_4;
fprintf('һ��ˮ�෵��:%f',h);
fprintf('һ��ˮ�������%f',V);
fprintf('һ��ˮ�ཬ�ܶ�:%f',roucm);
fprintf('һ����ˮ������kg��%f',W_c*1000);
fprintf('һ����ˮ����m^3:%f',Vw);
fprintf('һ��ʯӢɰ����kg��%f',W_c*1000/0.59*0.26);
fprintf('һ��΢������kg��%f',W_c*1000/0.59*0.15);
fprintf('ˮ�ұȣ�%f',m);
%fprintf('ˮ��������ƽ����/n');
%fprintf('�׸�ֱ��mm��%f',);
%fprintf('����/min��%f',);
%fprintf('�����L/s��%f',);




%%%%%
%����%
%����%
%%%%%

%�������þ������׹��ھ����׹��⾶,����
D_h=311.1*10^-3;
d=224.4*10^-3;
D=244.5*10^-3;
xx = input('���Enter a number: ');
yy = input('ѧ�ź���λ1234Enter a number: ');
h = 1500+(xx-1)*yy*3;

%ˮ����������K1:��������ϵ����rubber_stopper��pocket��
%
K1=1.1;
RS= 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huoˮ��
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


%ˮ����������
%��������
fai_600=56;fai_300=40;
n=3.32*log10(fai_600-fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);


%�����Ϊ3��
%��ŵ���������ٽ���ŵ����m/s;L/s
Re_c = 2100;
v_ci=(Re_c*800^n/(8000*roucm*d^n))^(1/(2-n));
v_ca=(0.2625*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
Q_ci=pi/40*d^2*v_ci;
Q_ca=pi/40*(D_h^2-D^2)*v_ca;
%�����ٽ������ٽ�����
Re_cp = 100;
v_ci_piston = (Re_cp*800^n*K/(8000*roucm*d^n))^(1/(2-n));
v_ca_piston = (0.0125*800^n*K/(roucm/(D_h^2-D^2)))^(1/(2-n));
Q_ca_piston = pi/40*(D_h^2-D^2)*v_ca_piston;
%�׹�Ь�ز�����ѹ��rouc:�꾮Һ�ܶ�h_cmˮ�෵��;h_�׹ܾ��pf=MPa;
h_cm=0;
pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);
%����ʩ����߱�ѹ
p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
p_ha=10^-3*9.81*(roudf*(h-h_cm)+roucm*h_cm);
%��Һ��ѹ��
delta_p = p_ha - p_hi;

x = input('����1Enter a number: ');

if x
    fi=0.03164/Re_c;
else
    a = (log10(n)+2.5)/50;
    b = (1.4-log10(n))/7;
    fi = a/Re^b;
end

%��������
p_fi = 2*h*roudf*v_ci^2*fi/d;
%��������
p_fa = 2*h*roudf*v_ca^2*fi/(D_h-D);

%summarise
p_max = p_ha +p_hi +p_fi +p_fa;


%%%עˮ��ʱ��
y = input('�ó�1234Enter a number: ');

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
z = input('�ó�����Enter a number: ');
V_m = (pi/4)*d^2*h;
T_1 = V/(Q_pump*z);
T_2 = 3;
T_3 = V_m/Q_pump;
T_4 = 3;
T = T_1+ T_2+ T_3+ T_4;
fprintf('����ˮ�෵��:%f',h);
fprintf('����ˮ�������%f',V);
fprintf('����ˮ�ཬ�ܶ�:%f',roucm);
fprintf('������ˮ������kg��%f',W_c*1000);
fprintf('������ˮ����m^3:%f',Vw);
fprintf('����ʯӢɰ����kg��%f',W_c*1000/0.59*0.26);
fprintf('����΢������kg��%f',W_c*1000/0.59*0.15);
fprintf('ˮ�ұȣ�%f',m);



%%%%%
%����%
%����%
%%%%%

%�������þ������׹��ھ����׹��⾶,����
D_h=215.9*10^-3;
d=121.36*10^-3;
D=139.7*10^-3;
h = 2021;

%ˮ����������K1:��������ϵ����rubber_stopper��pocket��
%
K1=1.1;
RS= 4;
pocket=4;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huoˮ��
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


%ˮ����������
%��������
fai_600=56;fai_300=40;
n=3.32*log10(fai_600-fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);


%�����Ϊ3��
%��ŵ���������ٽ���ŵ����m/s;L/s
Re_c = 2100;
v_ci=(Re_c*800^n/(8000*roucm*d^n))^(1/(2-n));
v_ca=(0.2625*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
Q_ci=pi/40*d^2*v_ci;
Q_ca=pi/40*(D_h^2-D^2)*v_ca;
%�����ٽ������ٽ�����
Re_cp = 100;
v_ci_piston = (Re_cp*800^n*K/(8000*roucm*d^n))^(1/(2-n));
v_ca_piston = (0.0125*800^n*K/(roucm/(D_h^2-D^2)))^(1/(2-n));
Q_ca_piston = pi/40*(D_h^2-D^2)*v_ca_piston;
%�׹�Ь�ز�����ѹ��rouc:�꾮Һ�ܶ�h_cmˮ�෵��;h_�׹ܾ��pf=MPa;
h_cm=0;
pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);
%����ʩ����߱�ѹ
p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
p_ha=10^-3*9.81*(roudf*(h-h_cm)+roucm*h_cm);
%��Һ��ѹ��
delta_p = p_ha - p_hi;

x = input('����1Enter a number: ');

if x
    fi=0.03164/Re_c;
else
    a = (log10(n)+2.5)/50;
    b = (1.4-log10(n))/7;
    fi = a/Re^b;
end

%��������
p_fi = 2*h*roudf*v_ci^2*fi/d;
%��������
p_fa = 2*h*roudf*v_ca^2*fi/(D_h-D);

%summarise
p_max = p_ha +p_hi +p_fi +p_fa;


%%%עˮ��ʱ��
y = input('�ó�1234Enter a number: ');

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
z = input('�ó�����Enter a number: ');
V_m = (pi/4)*d^2*h;
T_1 = V/(Q_pump*z);
T_2 = 3;
T_3 = V_m/Q_pump;
T_4 = 3;
T = T_1+ T_2+ T_3+ T_4;
fprintf('����ˮ�෵��:%f',h);
fprintf('����ˮ�������%f',V);
fprintf('����ˮ�ཬ�ܶ�:%f',roucm);
fprintf('������ˮ������kg��%f',W_c*1000);
fprintf('������ˮ����m^3:%f',Vw);
fprintf('����ʯӢɰ����kg��%f',W_c*1000/0.59*0.26);
fprintf('����΢������kg��%f',W_c*1000/0.59*0.15);
fprintf('ˮ�ұȣ�%f',m);





