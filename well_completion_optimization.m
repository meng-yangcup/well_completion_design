%ahead
%�������þ������׹��ھ����׹��⾶,���λm
D_h=444.5*10^-3;
d=313.6*10^-3;
D=339.7*10^-3;
h=260;
%ˮ����������K1:��������ϵ���� RS=rubber_stopper��pocket��
%
K1=1.1;
RS = 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huoˮ��
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

%ˮ����������
%��������
fai_600=56;fai_300=40;
n=3.32*log10(fai_600/fai_300);
K=0.511*fai_300/(511^n);
miu_p = fai_600-fai_300;
tao_0 = 0.511*(fai_300-miu_p);
miudf_p = 4.6;


%%%ѡ����������
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
        k = input('�ó�0123456Enter a number: ');
        if k==0
            continue
        else
            for j=i:3
                c1(j,:)=c(k,:);
            end
        end
    end
    p_pump = (sum(c1(:,4)));%MPa%Ϊ�򻯼�����Ϊ��ѹ�����������Ե���
    Q_pump = (sum(c1(:,3)))/1000*60;%m^3/min�����ЧΪ100% %
    V_m = (pi/4)*d^2*h;% volume of mud

    T_1 = V/(Q_pump);
    T_2 = 3;
    T_3 = V_m/Q_pump;
    T_4 = 3;
    T = T_1+ T_2+ T_3+ T_4;

%�����Ϊ2��������������
%%������ŵ����
%��������
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*(d*10)^n *vi^(2-n)/(800^n*K);

%������
    va = Q_pump*60/(pi/4*(D_h^2-D^2));
    Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
    %��ŵ���������ٽ���ŵ����m/s;L/s
        fprintf('ˮ������������ѡ������');
        fprintf('������Ϊ����');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
%    v_ci=(Re_i*800^n/(8000*roudf*d^n))^(1/(2-n));%����������
%    v_ca=(Re_a*800^n*K/(8000*roucm/(D_h-D)^n))^(1/(2-n));%��������
%    Q_ci=pi/40*d^2*v_ci;                         %��������
%    Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a < 100
        fprintf('������Ϊ����');
        fprintf('ˮ������������ѡ������');
        v_ci = Q_pump/(pi/4*d^2);
        v_ca = Q_pump/(pi/4*(D_h^2-D^2));
    %v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
    %v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
    %Q_ci=pi/40*d^2*v_ci;
    %Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('ˮ������������������ѡ��');
        continue;
    end

%�׹�Ь�ز�����ѹ��rouc:�꾮Һ�ܶ�h_cmˮ�෵��;h_�׹ܾ��pf=MPa;
    h_cm=0;
%pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:�����ܶ�
    rouf=100.0;
    pf=rouf* 9.81* (h-h_cm)*0.001;%MPa
%����ʩ����߱�ѹ
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%��Һ��ѹ��
    delta_p = p_ha - p_hi;
    if Re_a > 2100
        x = input('����1else����Enter a number: ');
        if x
            fi=0.03164/Re_a;
        else
            a = (log10(n)+2.5)/50;
            b = (1.4-log10(n))/7;
            fi = a/Re_a^b;
        end

%��������
        p_fi = 2*h*roudf*1000*v_ci^2*fi/d;
%��������
        p_fa = 2*h*roucm*1000*v_ca^2*fi/(D_h-D);
    else
        fi = 16/Re_i;
        fa = 24/Re_a;
    %%�����������ڻ�����������ܸ��ӣ���ô�򻯿��ܲ���
%��������
        p_fi = 2*h*roudf*1000*v_ci^2*fi/d;
%��������
        p_fa = 2*h*roucm*1000*v_ca^2*fi/(D_h-D);
    end
%summarise
    p_max = delta_p +p_fi +p_fa;
    if p_pump - p_max < pf
        fprintf('�ó�ѹ��ѡ�����');
        break;
    else
        fprintf('ѹ������������,������ѡ��');
    end
    
end

fprintf('һ��ˮ�෵��:%f\n',h);
fprintf('һ��ˮ�������%f\n',V);
fprintf('һ��ˮ�ཬ�ܶ�:%f\n',roucm);
fprintf('һ����ˮ������kg��%f\n',W_c*1000);
fprintf('һ����ˮ����m^3:%f\n',Vw);
fprintf('һ��ʯӢɰ����kg��%f\n',W_c*1000/0.59*0.26);
fprintf('һ��΢������kg��%f\n',W_c*1000/0.59*0.15);
fprintf('ˮ�ұȣ�%f\n',m);
%fprintf('ˮ��������ƽ����/n');
for i = 1:3
    fprintf('�׸�ֱ��mm��%f',c1(i:1));
    fprintf('����/min��%f',c1(i,2));
    fprintf('�����L/s��%f',c1(i,3));
    fprintf('���ѹMPa��%f',c1(i,4));
    fprintf('\n');
end
fprintf('עˮ��ʩ��ʱ�䣺%f\n',T);
fprintf('���꾮Һ�����%f\n',V_m);
fprintf('������%f\n',Q_pump);
fprintf('���ڷ��٣�%f\n',vi);
fprintf('���շ��٣�%f\n',va);

%%%%%
%����%
%����%
%%%%%


%�������þ������׹��ھ����׹��⾶,������¸�ֵ
D_h=311.1*10^-3;
d=224.4*10^-3;
D=244.5*10^-3;
xx = input('���Enter a number: ');
yy = input('ѧ�ź���λ1234Enter a number: ');
h = 1500+(xx-1)*50+yy*3;

%ˮ����������K1:��������ϵ����rubber_stopper�� pocket��
%
K1=1.1;
RS= 4;
pocket=0;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huoˮ��
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


%ˮ����������
%��������
%fai_600=56;fai_300=40;
%n=3.32*log10(fai_600-fai_300);
%K=0.511*fai_300/(511^n);
%miu_p = fai_600-fai_300;
%tao_0 = 0.511*(fai_300-miu_p);

%%%ѡ����������
pf=0;p_pump=1;p_max=0;
while pf<p_pump-p_max
    c1=zeros(3,4);
    for i=1:3
        k = input('�ó�0123456Enter a number: ');
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

%�����Ϊ2��������������
%%������ŵ����
%��������
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*d^n *vi^(2-n)/(800^n*K);
%������
    va = Q_pump/60/(pi/4*(D_h^2-D^2));
    Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
    %��ŵ���������ٽ���ŵ����m/s;L/s
        fprintf('������Ϊ����');
        fprintf('ˮ��������������ѡ������');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
        Q_ci=pi/40*d^2*v_ci;
        Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a < 100
        fprintf('������Ϊ����');
        fprintf('ˮ��������������ѡ������');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
        Q_ci=pi/40*d^2*v_ci;
        Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('��ǰ������ŵ��%f ',Re_a);
        fprintf('ˮ������������������ѡ��');
        continue;
    end

%�׹�Ь�ز�����ѹ��rouc:�꾮Һ�ܶ�h_cmˮ�෵��;h_�׹ܾ��pf=MPa;
    h_cm=0;
%pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:�����ܶ�
    rouf=10.0;
    pf=rouf* 9.81* (h-h_cm)/1000;%MPa
%����ʩ����߱�ѹ
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%��Һ��ѹ��
    delta_p = p_ha - p_hi;
    if Re_a > 2100

        x = input('����1else����Enter a number: ');
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
    %%�����������ڻ�����������ܸ��ӣ���ô�򻯿��ܲ���
    end
%��������
    p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;
%��������
    p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;
%summarise
    p_max = delta_p +p_fi +p_fa;
    if p_pump-p_max < pf
        fprintf('ѹ����������');
        
    else
        fprintf('ѹ��������������ѡ����');
    end
end
fprintf('����ˮ�෵��:%f\n',h);
fprintf('����ˮ�������%f\n',V);
fprintf('����ˮ�ཬ�ܶ�:%f\n',roucm);
fprintf('������ˮ������kg��%f\n',W_c*1000);
fprintf('������ˮ����m^3:%f\n',Vw);
fprintf('����ʯӢɰ����kg��%f\n',W_c*1000/0.59*0.26);
fprintf('����΢������kg��%f\n',W_c*1000/0.59*0.15);
fprintf('ˮ�ұȣ�%f\n',m);

%fprintf('ˮ��������ƽ����/n');
for i = 1:3
    fprintf('�׸�ֱ��mm��%f',c1(i:1));
    fprintf('����/min��%f',c1(i,2));
    fprintf('�����L/s��%f',c1(i,3));
    fprintf('���ѹMPa��%f',c1(i,4));
    fprintf('\n');
end
fprintf('עˮ��ʩ��ʱ�䣺%f\n',T);
fprintf('���꾮Һ�����%f\n',V_m);
fprintf('������%f\n',Q_pump);
fprintf('���ڷ��٣�%f\n',vi);
fprintf('���շ��٣�%f\n',va);

%%%%%
%����%
%����%
%%%%%

%�������þ������׹��ھ����׹��⾶,����
D_h=215.9*10^-3;
d=121.36*10^-3;
D=139.7*10^-3;
h = 2021;

%ˮ����������K1:��������ϵ����rubber_stopper�� pocket��
%
K1=1.1;
RS= 4;
pocket=4;
V = (pi/4)*(K1*h*(D_h^2-D^2)+d^2*RS+D^2*pocket);
%huoˮ��
rou_drycement=0.26*2.62+0.15*0.7+3.16*(1-0.15-0.26);
rou_drillingfluid = 1.18;
roudf=rou_drillingfluid;
%
K2=1.05;
roucm=1.82;rouw=1; m=0.515;Q=1.2;
q=roucm*rouw/(rouw+m*roucm);
W_c = K2* V *q;
Vw= m*W_c/rouw;

%%%ѡ����������
pf=0.1;p_pump=1;p_max=0;
while pf<p_pump-p_max
    c1=zeros(3,4);
    for i=1:3
        k = input('�ó�0123456Enter a number: ');
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

%�����Ϊ2��������������
%%������ŵ����
%��������
    vi = Q_pump/60/(pi/4*d^2); %m/s
    Re_i = 8000*roudf*d^n *vi^(2-n)/(800^n*K);

%������
    va = Q_pump/60/(pi/4*(D_h^2-D^2));
    Re_a = 8000*roucm*(D_h-D)^n *va^(2-n)/(800^n*K);
    if Re_a > 2100
    %��ŵ���������ٽ���ŵ����m/s;L/s
        fprintf('������Ϊ����');
        fprintf('����ѡ������');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
 %   Q_ci=pi/40*d^2*v_ci;
 %   Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    elseif Re_a < 100
        fprintf('������Ϊ����');
        fprintf('����ѡ������');
        v_ci=(Re_i*800^n/(8000*roucm*d^n))^(1/(2-n));
        v_ca=(Re_a*800^n*K/(roucm/(D_h-D)^n))^(1/(2-n));
%    Q_ci=pi/40*d^2*v_ci;
%    Q_ca=pi/40*(D_h^2-D^2)*v_ca;
    
    else
        fprintf('��ǰ��ŵ��Ϊ��%f',Re_a);
        fprintf('ˮ������������������ѡ��');
        continue;
    end

%�׹�Ь�ز�����ѹ��rouc:�꾮Һ�ܶ�h_cmˮ�෵��;h_�׹ܾ��pf=MPa;
    h_cm=0;
    %pf=10^-3*9.81*(h*roucm+ (h-h_cm)*roudf);rouf:�����ܶ�
    rouf=2.0;
    pf=rouf* 9.81* (h-h_cm);
%����ʩ����߱�ѹ
    p_hi=10^-3*9.81*(roudf*(h-RS)+roucm*RS);
    p_ha=10^-3*9.81*(roucm*(h-h_cm)+roucm*h_cm);
%��Һ��ѹ��
    delta_p = p_ha - p_hi;
    if Re_a > 2100

        x = input('����1else����Enter a number: ');
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
%%�����������ڻ�����������ܸ��ӣ���ô�򻯿��ܲ���
%��������
    p_fi = 2*h*roudf*1000*vi^2*fi/d*10^-6;
%��������
    p_fa = 2*h*roucm*1000*va^2*fi/(D_h-D)*10^-6;
    end
%summarise
    p_max = delta_p +p_fi +p_fa;

    if p_max - pf<pf
        fprintf('ѹ����������');
        break;
    else
        fprintf('ѹ��������������ѡ����');
    end
end
fprintf('����ˮ�෵��:%f\n',h);
fprintf('����ˮ�������%f\n',V);
fprintf('����ˮ�ཬ�ܶ�:%f\n',roucm);
fprintf('������ˮ������kg��%f\n',W_c*1000);
fprintf('������ˮ����m^3:%f\n',Vw);
fprintf('����ʯӢɰ����kg��%f\n',W_c*1000/0.59*0.26);
fprintf('����΢������kg��%f\n',W_c*1000/0.59*0.15);
fprintf('ˮ�ұȣ�%f\n',m);
%fprintf('ˮ��������ƽ����/n');
for i = 1:3
    fprintf('�׸�ֱ��mm��%f',c1(i:1));
    fprintf('����/min��%f',c1(i,2));
    fprintf('�����L/s��%f',c1(i,3));
    fprintf('���ѹMPa��%f',c1(i,4));
    fprintf('\n');
end
fprintf('עˮ��ʩ��ʱ�䣺%f\n',T);
fprintf('���꾮Һ�����%f\n',V_m);
fprintf('������%f\n',Q_pump);
fprintf('���ڷ��٣�%f\n',vi);
fprintf('���շ��٣�%f\n',va);





