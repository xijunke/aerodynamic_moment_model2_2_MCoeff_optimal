% function obj_function=objfunction_aeroM_stdev(x)
%% ��������: Solution of the aerodynamic moment
% �޸�ʱ�䡪��2014��12��23��,21:58����ת����ƫ�������е��ƫ��������,��Ťת������ƫ��C_maxy֮��
% ���ǳ�������ء�����Եڶ��ֽ��ٶ������
% 2014��6��20��,14:54:46
% 2014��6��20��,2:02:08
% ע�ⵥλ����(N.m)��*10^6��(mN.mm)=(uN.m)
% clear all;clc;
tic                             % Elapsed time is 99.240238 seconds.
x=[1, 1, 1, 1, 1.11, 1, 1, 1, 1, 1, 1, 1, 1, 1]; 
k_xaero=x(1);
k_xRot=x(2);
C_xRD=x(3);
k_am_x=x(4);
k_hinge_x=x(5);
k_inert_x=x(6);
k_inert_y=x(7);
k_inert_z=x(8);
k_ztrans=x(9);
k_zrot=x(10);
k_za=x(11);
C_zRD=x(12);
C_Ty1=x(13);
C_Ty2=x(14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ú�����ò��������������ϵ���ĺ���
% wing_para_output=zeros(1,16);
% wing_para_output=[F_ndTrans,Coeff_liftdragF_N,M_xaercoeff,I1y,Z_rnd, M_xrdcoeff,...
%     F_ndRot,F_yrotcoeff,M_xRotcoeff, I2y,...
%     I_xzam,I_xxam, I5z,I6z,I7y,M_zrdcoeff];
wing_para=wing_shape_fruitfly_sixteen_good();   %���ú���wing_shape_fruitfly;  % size(wing_para)
% (1) ƽ�����������������ز���
% F_ndTrans=wing_para(1,1);             % F_ndTrans=0.46391;  ������, ���ٻ���λΪmm^4
% Coeff_liftdragF_N=wing_para(1,2);  % Coeff_liftdragF_N=0.00682;  %��λ��mg*mm  %ƽ������������
M_xaercoeff=wing_para(1,3);              % M_xaercoeff=0.006038;   %��λ��: mg.mm^2   % ƽ�������������ز��������Ƴ�ƽ���µ�չ����
% C_aver1=M_xaercoeff/Coeff_liftdragF_N;  % C_aver1=0.8854;
I1z=wing_para(1,4);                             % I1y=0.016158   % ��λ�� mg.mm^2            % ƽ�������������ز��������Ƴ�ƽ���µ�������
% Z_rnd=wing_para(1,5);                     % Z_rnd=0.16265;  ������, ���ٻ���λΪmm
M_xrdcoeff=wing_para(1,6);                % M_xrdcoeff=0.0001839; % ��λ��mg.mm^2 %ת�������������ز������Ƴ�ƽ���µ�չ����
% (2) ת�����������������ز���
% F_ndRot=wing_para(1,7);                % F_ndRot=0.74847;  ������, ���ٻ���λΪmm^4
% F_yrotcoeff=wing_para(1,8);               % F_yrotcoeff =0.003243;  % ��λ�� mg.mm   % ת������������
M_xRotcoeff=wing_para(1,9);             % M_xRotcoeff=0.002871;   % ��λ�� mg.mm^2 % ת��������������ϵ�������Ƴ�ƽ���µ�չ����
% C_aver2=M_xRotcoeff/F_yrotcoeff    % C_aver2=0.8854;
I2z=wing_para(1,10);                          % I2y=0.006943;        % ��λ�� mg.mm^2            % ת�������������ز��������Ƴ�ƽ���µ�������
% (3) �����������������ز���
I_xzam=wing_para(1,11);                    % I_xzam =0.001424  % ��λ�� mg.mm^2  % �������������ز��������Ƴ�ƽ���µ�չ����
I_xxam=wing_para(1,12);                    % I_xxam =0.000338  % ��λ�� mg.mm^2  % �������������ز��������Ƴ�ƽ���µ�չ����
% I5y=wing_para(1,13);                          % I5z=0.0050926   % ��λ�� mg.mm       % ������������������������
% I6y=wing_para(1,14);                          % I6z=0.00077164    % ��λ�� mg.mm       % ������������������������
I7z=wing_para(1,15);                          % I7y=0.0109056;      % ��λ�� mg.mm^2   % �������������ز��������Ƴ�ƽ���µ�������
M_zrdcoeff=wing_para(1,16);             % M_zrdcoeff=0.001169; % ��λ�� mg.mm^2 % ת�������������ز������Ƴ�ƽ���µ�������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I5y=wing_para(1,13);                       % I5z=0.0050945   % ��λ�� mg.mm       % ������I3yӦ�ø�ΪI5y
% I6y=wing_para(1,14);                       % I6z=0.0011         % ��λ�� mg.mm       % ������I4yӦ�ø�ΪI6y
% I7z=wing_para(1,15);                       % I7y=0.0109;        % ��λ�� mg.mm^2   % ������I5zӦ�ø�ΪI7z
% I1z=wing_para(1,4);                         % I1y=0.0162         % ��λ�� mg.mm^2    % ������I7zӦ�ø�ΪI1z
% I2z=wing_para(1,10);                       % I2y=0.0069;         % ��λ�� mg.mm^2    % ������I6zӦ�ø�ΪI2z    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���ó���˶�ѧ-���˶����ɺͼ��ι���(AOA)������
% wing_m_output=[t',phi',psi',alpha',dphi',dpsi',ddphi',ddpsi',C_L',C_D',C_N1',C_T'];
wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();      %���ú���kenimatics_wing_and_AoA;  % size(wing_kenimatics)  % (1000,12)
t=wing_kenimatics(:,1);               % ��λ��ms
phi=wing_kenimatics(:,2);           % �Ĵ�ǡ�����λ��rad
psi=wing_kenimatics(:,3);            % �Ĵ�ǡ�����λ��rad
% alpha=wing_kenimatics(:,4);       % alpha=atan2(omega_z,-omega_y);  % ���ι��ǡ���������   %������������и�
dphi=wing_kenimatics(:,5);         % ��λ��rad/s
dpsi=wing_kenimatics(:,6);         % ��λ��rad/s
ddphi=wing_kenimatics(:,7);      % ��λ��rad/s^2
ddpsi=wing_kenimatics(:,8);       % ��λ��rad/s^2
% C_L=wing_kenimatics(:,9);          
% C_D=wing_kenimatics(:,10);     
C_N1=wing_kenimatics(:,11);   
C_T=wing_kenimatics(:,12); 
C_N=C_N1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڳ�����ϵ����������������������ϵ��: Omega1; Omega2; Omega3; Omega4;
% ���ڵ�2�ַ�������������ء�������ϵ��(����ϵԭ��Ϊ�����O'): ��������ϵ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413; %ǰԵ����Ϻ���
% % yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282; %��Ե����Ϻ���
% Omega1= 0.01615796930;           %��λ�� mg.mm^2  ���������ܻ��������Ҫ�˶Լ��㡪��ƫС
% Omega2 = 0.0004196574748;      %��λ�� mg.mm^2  ���������ܻ��������Ҫ�˶Լ���
% Omega3=0.006942725660;          %��λ�� mg.mm^2   ���������ܻ��������Ҫ�˶Լ���
% Omega4=0.007031892405;          %��λ�� mg.mm^2   ���������ܻ��������Ҫ�˶Լ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������ϵ�µĽ����ʺͽǼ����ʡ����������������������Գ�2DOF�˶�
% f=188.7;    % Hz
% T=1/f;        
% �������ϵ�µĽ��ٶ�
omega_x=dpsi;                      % չ��
omega_y=dphi.*sin(psi);        % ����(��ʼ����)
omega_z=dphi.*cos(psi);       % ������
omega_h=dphi;       % �����Ľ��ٶ�% omega_h=-sqrt(omega_y^2+omega_z^2);  % ���ַ���˳ʱ��
% �������ϵ�µĽǼ��ٶȡ����������������ļ���
domega_x=ddpsi;
domega_y=ddphi.*sin(psi)+dphi.*dpsi.*cos(psi);
domega_z=ddphi.*cos(psi)-dphi.*dpsi.*sin(psi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������Ǻ͵��������ٶȵļ��㡪��ȡ�����ٶ�����ڸ�������ٶȣ���������Ӹ��ţ�
v_y_nonr=-omega_z;     % v_y=r*dphi*cos(psi)
v_z_nonr=omega_y;       % v_z=-r*dphi*sin(psi)
alpha2=atan2(-v_y_nonr,v_z_nonr);   % ��ȷ����ע�������ĵ�alpha2=atan2(omega_z,-omega_y)*180/pi; ��ͬ
% % ����alpha2=atan2(cot(psi))=atan2(cot(pi/2-alpha))=atan2(tan(alpha)); %����atan2������������ֵ��������alpha>pi/2ʱ
V_nonr=sqrt(v_y_nonr.^2+v_z_nonr.^2); % ���������ٶ�V_nonr=omega_h=dphi;   % ��λ�� rad/s
% figure(11)
% plot(t/T,alpha*180/pi,'r-',t/T,alpha2*180/pi,'b-.','LineWidth',2)  
% xlabel('\itNormalized time')
% ylabel('\it\alpha & AoA  (deg)')
% legend('\alpha(t)','AoA(t)')
% title('����AoA��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һģ�顪��չ��Ťת�ᡪ���������ط���%%%%
%%%%���߲���, �ֱ��ǣ�%%%%%%%%%%%%%%��λ: (mN.mm)��(uN.m)
%%%%��һ���֡���ƽ��������������������%%%%%%
%%%%�ڶ����֡���ת��������������������%%%%%%
%%%%�������֡���ת��������������%%%%%%%%%
%%%%���Ĳ��֡���ת������������%%%%%%%%%%
%%%%���岿�֡���Ťת�����ĵ��Իظ�����%%%%%%
%%%%�������֡������������%%%%%%%%%%%%
%%%%���߲��֡�����������������%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���֡���ƽ�������������������ء�Ťת�ᡪ��λ: (mN.mm)��(uN.m)
% �������Ϊ��5��
% C_avereff;   R_wingeff;   F_nd;   Y_rcpnd;         �����������������������Գ���ò������
% alpha2;   omega_h;   C_N(alpha2);                      �����������������������Գ�2DOF�˶�����,����������
% ƽ��Ťת�������ء�����ת����������
% M_xtrans=(-sign(alpha2).*Rou.*omega_h.^2.*C_N.*C_avereff^2*R_wingeff^3*F_nd.*Y_rcpnd/2)*10^(-15); %ԭ�Ĺ�ʽN.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % (1) ����1����ƽ������������ת�������ء���ѹ��λ��Ťת��֮��
% % F_ytran=-sign(alpha2).*abs(C_N).*V_nonr.^2*Coeff_liftdragF_N*10^(-3);   % ��λ��rad^2*s^-2*mg*mm=10^(-9)N=10^(-9)*10^6uN
% % ����ĵ�λ�� (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
% Y_rcpnd=COP_Ycpnd2_fruitfly(alpha2);  % �����ת���������ء������ú���COP_Ycpnd2_fruitfly��⾻ѹ�ĵ�������λ��Y_rcpnd
% % Y_rcpndaver=mean(Y_rcpnd)    % Y_rcpndaver=0.1298(old); % Y_rcpndaver= - 0.0271;
% M_xtrans1=-sign(alpha2).*abs(C_N).*V_nonr.^2.*Y_rcpnd*M_xaercoeff*10^(-3);    % ��ת����������һ��ʼ����ʱ���
% Z_trans1=M_xtrans1./F_ytran;    % ������������
% % Z_transaver1=mean(Z_trans1)   % Z_transaver1=0.1149(old); % Z_transaver1= - 0.0240;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) ����2����ƽ������������ת�������ء���ѹ��λ��Ťת��֮�󡪡�Ťת�ᡪ�ڶ�������������صķ���
% ���������������Ϊ��10�� %% ��λ: (mN.mm)��(uN.m)
% Omega1; Omega2; Omega3; Omega4;  ��ǰ�ġ����������������������Գ���ò������
%%%%%%%%%%%%%%%%%%%%%%%%%
% % alpha0=pi/4;     % ���Ǽ���Ϊ45�㣬�����ת��������������ѹ��λ�á�����ʱ�������ù��Ǻ���������XXX
% d_cprnd=0.82*abs(alpha2)/pi+0.05;                  %����Ϊpi/4ʱ����ת���������ص�����ѹ��λ��;
% % M_rotx=((Omega3*d_cprnd-Omega4).*sign(dphi).*dphi.^2.*C_N)*(10^-12*10^6); % ��������Ťת��(x��)��������Pitching Motion
% % mg.mm^2*rad*s^-2=(10^-12)kg*m*s^-2*m=10^-12N*m=10^-3 uN*mm
% % k_xtranscirc=0.716;     % ��ϵ�������⡪��XXX
% k_xtranscirc=0.1;     % ��ϵ��������
% % uN*mm % ��������Ťת��(x��)�������׸�������ע������(Omega3*d_cprnd-Omega4)��0
% M_xtranscirc=-k_xtranscirc*sign(alpha2).*(Omega3*d_cprnd-Omega4).*dphi.^2.*C_N*(10^-3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(12)  
% plot(t/T,M_xtrans1,'b-',t/T,M_xtranscirc,'g-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{trans,x} & M_{transcirc,x} )  (uN.mm)')
% legend('M_{trans,x}','M_{transcirc,x}')
% title('չ��(x-axis): ƽ������������ת�����������ء����Աȷ���')
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) ����3���������ת���������ء�����ѹ��λ��
% k_xaero=1; 
N=length(t);
M_xtrans=zeros(N,1);
Y_rcpnd_trans=zeros(N,1);
for i=1:1:N
    Y_rcpnd_trans(i,1)=COP_Ycpnd2_TransCirc(alpha2(i,1));  % ���ú���COP_Ycpnd2_TransCirc��⾻ѹ�ĵ�������λ��Y_rcpnd; % ��������
    %  Y_rcpnd_trans=abs(Y_rcpnd_trans(i,1));  % ��
    % ����ĵ�λ�� (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
    %  M_xaercoeff=0.0060;   %��λ��: mg.mm^2
    % ƽ��������������ת����������һ��ʼ����ʱ���
     % M_xtrans(i,1)=-k_xaero*sign(alpha2(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans(i,1)*M_xaercoeff*10^(-3);      
     M_xtrans(i,1)=k_xaero*sign(alpha2(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans(i,1)*M_xaercoeff*10^(-3);  
end
% % M_xrotcirc(i,1)=k_xRot*C_R*omega_x(i,1).*abs(omega_h(i,1)).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3);
% % F_ytran=-sign(alpha2).*abs(C_N).*V_nonr.^2*Coeff_liftdragF_N*10^(-3); 
% Z_trans=M_xtrans./F_ytran;    % ������������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ƽ��������������������
% figure(13)  
% % plot(t/T,M_xtrans1,'b-',t/T,M_xtranscirc,'g-',t/T,M_xtrans,'r-','LineWidth',2)
% plot(t/T,M_xtrans1,'b-',t/T,M_xtrans,'g-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{trans1,x} & M_{transcirc,x} & M_{trans,x}  (uN.mm)')
% % legend('M_{trans1,x}','M_{transcirc,x}','M_{trans,x}')
% legend('M_{trans1,x}','M_{trans,x}')
% title('չ��(x-axis): ƽ������������ת�����������ء����Աȷ���')
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ����֡���ת�������������������ء���Ťת��
% k_xRot=1;  
C_R=1.55;    % ��ϵ�������޸�
% % C_R=2.58;
% % M_xtrans(i,1)=k_xaero*sign(alpha2(i,1)).*abs(C_N(i,1)).*V_nonr(i,1).^2.*Y_rcpnd_trans(i,1)*M_xaercoeff*10^(-3);  
% % F_yrot=C_R*omega_x.*abs(omega_h)*F_yrotcoeff*10^(-3); % ��λ��rad^2*s^-2*kg*m=10^(-9)N=10^(-9)*10^6uN
% F_yrot=C_R*omega_x.*abs(V_nonr)*F_yrotcoeff*10^(-3);      % ��λ��rad^2*s^-2*kg*m=10^(-9)N=10^(-9)*10^6uN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_xRotcoeff=k_xRot*C_R*M_xRotcoeff;
% M_xrotcirc1=omega_x.*abs(omega_h).*Y_rcpnd*M_xRotcoeff*10^(-3);
% Z_rot1=M_xrotcirc1./F_yrot;           % ������������
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(t);
M_xrotcirc=zeros(N,1);
Y_rcpnd_rot=zeros(N,1);
for i=1:1:N
    % ת��������Ťת���������ء������ú���COP_Ycpnd2_RotCirc��⾻ѹ�ĵ�������λ��Y_rcpnd_rot
    % ת�����������ġ�ѹ�ķֲ�����Dickinson���� or ѹ�������ҵ� or ѹ����c(r)/4����Ťת������
    Y_rcpnd_rot(i,1)=COP_Ycpnd2_RotCirc(alpha2(i,1));  % ѹ�ķֲ�����Dickinson����
   % Y_rcpnd_rot=abs(Y_rcpnd_rot(i,1));
   % ����ĵ�λ�� (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
   %  M_xrotcirc(i,1)=k_xRot*C_R*sign(alpha2(i,1)).*omega_x(i,1).*abs(V_nonr(i,1)).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3); % ת��������Ťת����������  
   M_xrotcirc(i,1)=-k_xRot*C_R*omega_x(i,1).*abs(omega_h(i,1)).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3);    % ת��������Ťת����������
end
% Z_rot=M_xrotcirc./F_yrot;           % ������������
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(14)  % ת���������ط���, ת���������ط�����ת��������Ťת����������
% hold on
% plot(t/T,M_xrotcirc1,'r-',t/T,M_xrotcirc,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{rotcirc,x})  (uN.mm)')
% legend('M_{rotcirc1,x}','M_{rotcirc,x}')
% title('չ��(x-axis): ת��������Ťת������������ʱ��ı仯����')
% grid on
% % axis([0.9,4.05,-inf,inf])
% % set(gca,'XTick',(1:0.1:4.05))
% axis([0.9,3,-inf,inf])
% set(gca,'XTick',(1:0.1:3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �����湥�Ǳ仯������ƽ��������ת��������������ѹ��
% % Z_rot=M_xrotcirc./F_yrot;           % ������������
% % Z_trans=M_xtrans./F_ytran;        % ������������
% Z_circ=(M_xtrans+M_xrotcirc)./(F_ytran+F_yrot);    % ������������
% % size(Y_rcpnd_trans)
% % size(Y_rcpnd_rot)
% % x_cop=[t,alpha2*180/pi,Y_rcpnd,Z_trans1,Y_rcpnd_trans,Z_trans,Y_rcpnd_rot,Z_rot,Z_circ]; 
% % xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\aerodynamic_moment\aerodynamic_moment_model2_2\x_cop.xlsx',x_cop,'sheet1','A1:I2000');
% % (a) ƽ�������������ص�����ѹ�ĵ�ƽ��ֵ
% % ������ú���COP_Ycpnd2_fruitfly��⾻ѹ�ĵ�������λ��Y_rcpnd
% Y_rcpndaver=mean(Y_rcpnd)                     % Y_rcpndaver=-0.0271;
% % ����Z_trans1=M_xtrans1./F_ytran;          % ������������
% Z_transaver1=mean(Z_trans1)                    % Z_transaver1=-0.0240;
% % ������ú���COP_Ycpnd2_fruitfly��⾻ѹ�ĵ�������λ��Y_rcpnd;
% Y_rcpnd_transaver=mean(Y_rcpnd_trans)   % Y_rcpnd_transaver=0.0144;%Y_rcpnd_transaver=trapz(t,Y_rcpnd_trans)/(3*T)
% % ����Z_trans=M_xtrans./F_ytran;               % ������������
% Z_transaver=mean(Z_trans)                         % Z_transaver= 0.0127;
% % (b) ת�������������ص�����ѹ�ĵ�ƽ��ֵ
% % ������ú���COP_Ycpnd2_RotCirc��⾻ѹ�ĵ�������λ��Y_rcpnd_rot
% Y_rcpnd_rotaver=mean(Y_rcpnd_rot)  % Y_rcpnd_rotaver=0.0133; % Y_rcpnd_rotaver=trapz(t,Y_rcpnd_rot)/(3*T) 
% % ����Z_rot=M_xrotcirc./F_yrot;                % ������������
% Z_rotaver=mean(Z_rot)                              % Z_rotaver=0.0118;
% % (c) �����������ص�����ѹ�ĵ�ƽ��ֵ
% % ����Z_circ=(M_xtrans+M_xrotcirc)./(F_ytran+F_yrot);    % ������������
% Z_circaver=mean(Z_circ)                             % Z_circaver=0.0092;
%% ƽ��������ת��������������������ѹ����ʱ��仯����Z_rot1
% figure(15) 
% plot(t/T,Y_rcpnd,'b:',t/T,Z_trans1,'g-.',t/T,Y_rcpnd_trans,'b-',t/T,Z_trans,'g-',t/T,Y_rcpnd_rot,'r:',t/T,Z_rot1,'r-.',t/T,Z_rot,'r-',t/T,Z_circ,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('ƽ��������ת��������������������ѹ����Ťת��֮��ľ���  (mm)')
% legend('Y_{rcpnd}(t)','Z_{trans,1}(t)','Y_{rcpnd,trans}(t)','Z_{trans}(t)','Y_{rcpnd,rot}(t)','Z_{rot1}(t)','Z_{rot}(t)','Z_{circ}(t)');
% title('ƽ��������ת������������������ѹ����ʱ��ı仯����')
% grid on
% % hold on
% % % plot(t/T,F_ytran,'r-',t/T,F_yrot,'g-',t/T,F_ytran+F_yrot,'k-',t/T,M_xtrans+M_xrotcirc,'k:','LineWidth',2)
% % plot(t/T,F_ytran+F_yrot,'k-',t/T,M_xtrans+M_xrotcirc,'k:','LineWidth',2)
% axis([0.9,4.05,-inf,inf])
% % axis([0.9,4.05,-0.3,0.3])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������֡���ת�������������ء�Ťת�ᡪ��λ: (mN.mm)��(uN.m)
% �����������ء����������Ϊ��2��
% C_avereff;   R_wingeff;    Z_rnd;     �����������������������Գ���ò������
% omega_x;                                      �����������������������Գ�2DOF�˶�������
% (1) ����1����ת��Ťת�������ء���ת��������������
% M_xrd=(-omega_x.*abs(omega_x)*Rou*C_RD*C_avereff^4*R_wingeff*Z_rnd/2)*10^(-15); %ԭ�Ĺ�ʽ N.m
% C_xRD=1;  
% C_RD=5.0;     %ת������������ϵ��:C_RD=C_rd��(3,6); ����C_rd=C_dmax=3.4;
% C_RD=2;  
% C_RD=1;  % �������������õ�ϵ��
% uN.mm   % (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
M_xrd=-C_xRD*omega_x.*abs(omega_x)*M_xrdcoeff*10^(-3);   % ��������һ��ʼ��˳ʱ��� % ע���������ȡ����omega_x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) ����2���� Ťת���ٶȲ����ġ���ת���������ء����ڶ�������������صķ���
% ���������������Ϊ��10�� %% ��λ: (mN.mm)��(uN.m)
% Omega1; Omega2; Omega3; Omega4;  ��ǰ�ġ����������������������Գ���ò������
% M_rdampx=-(Omega2*sign(dpsi).*dpsi.^2*C_RD)*(10^-12*10^6);  % ��������Ťת��(x��)��������Pitching Motion
% mg.mm^2*rad*s^-2=(10^-12)kg*m*s^-2*m=10^-12N*m=10^-3 uN*mm
% k_xrotdamp=0.3775;  % ��C_RD=5ʱ;  ��ϵ����������M_xrdƥ��
% k_xrotdamp=0.438;
% M_xrotdamp=-k_xrotdamp*sign(dpsi).*Omega2.*dpsi.^2*C_RD*(10^-3);  % uN*mm % ��������Ťת��(x��)���������׸���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(16)  
% hold on
% plot(t/T,M_xrd,'k-',t/T,M_xrotdamp,'rd','LineWidth',1.5)
% xlabel('\itNormalized time')
% ylabel('\itM_{rd,x}  & M_{rdamp,x} )  (uN.mm)')
% legend('M_{rd,x}','M_{rdamp,x}')
% title('չ��(x-axis): ת���������ء����Աȷ���') % ��������������ʱ��ı仯����
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �Աȷ�������ƽ�����������������غ�ת����������֮�͵Ĳ���
% figure(17)  % ���ּ��㷽����õı���Ťת���صĶԱ� 
% M_xpassiverot=M_xtrans1+M_xrd;   % ����Ťת�������غ�3�ֻ��Ƶ���������
% M_totalx=M_xtranscirc+M_xrotdamp;  % ��ƽ��: չ��Ťת��(x��)
% plot(t/T,M_xpassiverot,'r-',t/T,M_totalx,'b-','LineWidth',2)    
% xlabel('\itNormalized time')
% ylabel('\itM_{passiverot,x}  &  M_{total,x}  (mN.mm)')      
% legend('M_{passiverot,x}','M_{total,x}')
% title('����Ťת�������ص��ܺ�(չ��:x-axis)����������ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���Ĳ��֡���ת������������+ƽ�����������ء�Ťת�ᡪ��λ: (mN.mm)��(uN.m)
% ���������������������Ϊ��6�������������������ġ�ѹ�������ҵ㡪Ťת������
% I_xzam;   I_xxam;       �����������������������Գ���ò������
% omega_x;   omega_z;   domega_x;   domega_y;   �����������������������Գ�2DOF�˶������ʺͽǼ�����
% ������Ťת�������ء������ݹ�Ӭ����޸ĵĽ��ٶȺͽǼ��ٶ�
% M_xam=(-I_xyam*(domega_y-omega_x*omega_z)-I_xxam*domega_x)*10^(-12); %N.m-ԭʼ��ʽ����������ֵ����,��Ҫת����
% ����ĵ�λ��: mg.mm^2*(rad/s)^2=mg.mm/s^2.mm=10^(-3) uN.mm; %ע��������������Ϊ����������domega_x=ddpsi;  
% k_am_x=1;   % ����������ϵ��; 2.35�ǽ������ޣ� psi_max =32.9386
% F_yadd1=(-I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3); 
% ����������һ��ʼ��˳ʱ��ģ��ܿ���Ϊ��ʱ�룬�������ƽ�����ٺ�ת���������㶨���ٽ�����ƽ�����٣�ת������
% M_xam=k_am_x*(-I_xzam*(domega_z-omega_x.*omega_y)-I_xxam*domega_x)*10^(-3);  % uN.mm 
M_xam=k_am_x*(-I_xzam*(domega_z+omega_x.*omega_y)-I_xxam*domega_x)*10^(-3);   % ��ʼ��ʱ��(-)
% c_zcopnd_add=M_xam./F_yadd1;                              % ������������
% % c_zcopnd_addaver=trapz(t,c_zcopnd_add)/(3*T); % c_zcopnd_addaver =-0.3058;
% c_zcopnd_addaver=mean(c_zcopnd_add)                 % c_zcopnd_addaver=-0.3058; 
% figure(18)                         
% plot(t/T,M_xam,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{am,x}   (uN.mm)')
% legend('M_{am,x}')
% title('����������������ʱ��ı仯����')   % ����������������ʱ��ı仯����
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���岿�֡��� Ťת�����ظ����ء���Ťת��
% �������������Paremeter of wing hinge
L_h1=70e-006;     % L_h2=175e-006;   % length of wing hinge  \um�����Ѿ����㵽m
W_h=1.8e-003;                                     % width of wing hinge  \mm�����Ѿ����㵽m
t_h=7.6e-006;                                       % thickness of wing hinge \um�����Ѿ����㵽m
E_h=2.5e009;                                        % lmodulus of elasticity for wing hinge  \Gpa ע�ⵥλ
%����ĵ�λ�� N/m^2*m^3*m/m=N.m=kg*m/s^2.m=mg.mm/s^2.mm:[10^12]=10^6(mN.mm)=10^9(uN.mm)
% k_hinge_x=1.11;      % psi_max =70.4552 @ k_am=1;  % �����Ťת�Ƿ�ֵ��ʵ��Ťת�Ƿ�ֵ
k_h=k_hinge_x*E_h*t_h^3*W_h/(12*L_h1)*10^9;  % uN.mm;  % mg.mm/s^2.mm
M_hinge=-k_h*psi*10^(-3);    % k_hΪrotational stiffness of the passive hinge  % �����ظ����ؿ�ʼ��˳ʱ���
% figure(19)
% plot(t/T,M_hinge,'c-')  
% xlabel('\itNormalized time')
% ylabel('\itM_{hinge}   (uN.mm)')      
% legend('M_{hinge}')
% title('�����ظ�����(չ��:x-axis)��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������֡�����������ء���Ťת�ᡪ����ƽ�������µ����ط�������%������û��Ҫ����
m_wing=2.4*10^-9;                                             % kg
% g=9.821;                                                               % m/s^2=N/kg
xr=0.3289*10^-3;                                                %  \m       % x-root offset 
x_com=xr+1.920243385*10^-3;                         %  \m       % ���ĵ�չ������
% z_com=-0.149785466+0.636=0.486215*10^-3;% ��Ťת����������
z_com=0.149785466*10^-3;                               %  \m       % ��Ťת����������
% d_com=z_com;                                                    % \m        % ����������
% r_cog=x_com;
% M_weight_x=m_wing*g*d_com*sin(psi)*10^9;   % ��λ��: 10^9(uN.mm)  % �����ؿ�ʼ��˳ʱ���
% M_weight_y=m_wing*g*r_cog*cos(psi)*10^9;
% M_weight_z=-m_wing*g*r_cog*sin(psi)*10^9;
% figure(20) 
% plot(t/T,M_weight_x,'r-',t/T,M_weight_y,'b-',t/T,M_weight_z,'g-','LineWidth',2) 
% xlabel('\itNormalized time')
% ylabel('��������������ط���M_{inert,x}(t) & M_{inert,y}(t) & M_{inert,z}(t) (uN.mm)'); 
% legend('M_{weight,x}(t)','M_{weight,y}(t)','M_{weight,z}(t)');  
% title('��������������ط���M_{inert,x}(t),M_{inert,y}(t) & M_{inert,z}(t)��ʱ��ı仯')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���߲��֡����������������ء���Ťת��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_inert_y=-m_wing*(domega_z*x_com-domega_x*z_com+omega_y.*(omega_x*x_com+omega_z*z_com))*10^6;  % ������������������
F_inert_z=-m_wing*(-domega_y*x_com-omega_y.^2*z_com+omega_x.*(omega_z*x_com-omega_x*z_com))*10^6; % ������������������
% ��������������
M_inert_x=-k_inert_x*z_com*F_inert_y*10^3;  % չ�򡪡���������������;  ��ʼ��ʱ��(-)
M_inert_z=-k_inert_y*x_com*F_inert_y*10^3;  % ���򡪡���������������;  ��ʼ��ʱ��(-)
M_inert_y=k_inert_z*x_com*F_inert_z*10^3;   % ���򡪡���������������;  ��ʼ˳ʱ��(-)
% figure(21) 
% plot(t/T,M_inert_x,'k-',t/T,M_inert_z,'r-',t/T,M_inert_y,'g-','LineWidth',2) 
% xlabel('\itNormalized time')
% ylabel('�����������������ط���M_{inert,x}(t) & M_{inert,z}(t) & M_{inert,y}(t) (uN.mm)'); 
% legend('M_{inert,x}(t)','M_{inert,z}(t)','M_{inert,y}(t)');  
% title('�����������������ط���M_{inert,x}(t) & M_{inert,z}(t) & M_{inert,y}(t)��ʱ��ı仯')  
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ��1�ַ�����⡪������ϵ�£�����Ťת�������ص��ܺ� M_xpassive_rotation; ��������Ťת��(x��)
%������ϵ��: ����Ťת��������{��������Ťת��(x��)}: 
% % M_xpassiverot=M_xtrans1+M_xrd+M_xam;                % ����Ťת�������غ�3�ֻ��Ƶ���������
% M_xaero1=M_xtrans1+M_xrd+M_xam;                           % ����Ťת�������غ�3�ֻ��Ƶ���������
M_xaero=M_xtrans+M_xrd+M_xam+M_xrotcirc;            % ����Ťת�������غ�4�ֻ��Ƶ���������
M_xtotal=M_xaero+M_inert_x;
%%%%%%%%%%%%%%%%%%%%%%%%%
% M_xtotal1=M_xtrans1+M_xrd+M_xam-M_hinge;                                        % -M_inert_x-M_weight_x����������û��Ҫ����
% M_xtotal=M_xtrans+M_xrd+M_xam+M_xrotcirc-M_hinge;                          % -M_inert_x-M_weight_x����������û��Ҫ����
% M_xtotal2=M_xtrans+M_xrd+M_xam+M_xrotcirc-M_hinge+M_inert_x;           % -M_weight_x����6��
% M_xtotal1=M_xtrans1+M_xrd+M_xam-M_hinge+M_inert_x;                        % -M_weight_x����5���������û��Ҫ����
% M_xtotal1=M_xtrans1+M_xrd+M_xam-M_hinge+M_inert_x+M_weight_x;                    %  6���������û��Ҫ����
% M_totalx=M_xtrans+M_xrd+M_xam+M_xrotcirc-M_hinge+M_inert_x+M_weight_x;     %  7���������û��Ҫ����
%% (1)�����������ѹ�ļ���ƽ����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(22)  % ƽ��������������, ת���������ء�ת�������������غ�����������-��Ťת����������
% plot(t/T,M_xtrans1,'r-',t/T,M_xrd,'g-',t/T,M_xam,'b-',t/T,M_xaero1,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{aero1,x}  ( including M_{trans1,x}  & M_{rd,x}  &  M_{am,x} &)  (uN.mm)')
% legend('M_{trans1,x}','M_{rd,x}','M_{am,x}','M_{aero1,x}')
% title('չ��(x-axis)-ƽ��������������,ת����������,���������غ�����������-��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
% % hold on
% % plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
%% (2)����Ƭ������ѹ�ļ���ƽ����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(23)  % ƽ��������������, ת���������ء�ת�������������غ�����������-��Ťת����������
% plot(t/T,M_xtrans,'r-',t/T,M_xrd,'g-',t/T,M_xrotcirc,'m-',t/T,M_xam,'b-',t/T,M_xaero,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{aero,x}  ( including M_{trans,x}  & M_{rd,x}  & M_{rotcirc,x} & M_{am,x})  (uN.mm)')
% legend('M_{trans,x}','M_{rd,x}','M_{rotcirc,x}','M_{am,x}','M_{aero,x}')
% title('չ��(x-axis)-ƽ��������������,ת����������,ת��������������,���������ؼ�����������-��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
% % hold on
% % plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(24)  % չ��(x-axis)����������ʱ��ı仯����
% plot(t/T,M_hinge,'c-',t/T,M_inert_x,'y',t/T,M_xaero1,'m-',t/T,M_xtotal1,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{total1,x}  ( including M_{x,aero1} & M_{hinge} & M_{inert,x})  (uN.mm)')
% legend('M_{hinge}','M_{inert,x}','M_{aero1,x}','M_{total1,x}')
% title('չ��(x-axis)-��������,�����ظ����غͳ��������������ؼ�������-��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
% hold on
% set(gca,'XTick',(0.9:0.1:4.05))
% % hold on
% % plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(25)  % չ��(x-axis)����������ʱ��ı仯����
% plot(t/T,M_hinge,'c-',t/T,M_inert_x,'y',t/T,M_xaero,'m-',t/T,M_xtotal,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{total,x}  ( including M_{x,aero} & M_{hinge} & M_{inert,x})  (uN.mm)')
% legend('M_{hinge}','M_{inert,x}','M_{aero,x}','M_{total,x}')
% title('չ��(x-axis)-��������,�����ظ����غͳ��������������ؼ�������-��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
% hold on
% set(gca,'XTick',(0.9:0.1:4.05))
% % hold on
% % plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure(26)        % ͼ1�����Ĵ��,Ťת�Ǻͼ��ι���AOA
% plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
% % xlabel('\itNormalized time')
% % ylabel('\itAngle (��)')
% % legend('\it\phi(t)','\it\psi(t)','\it\alpha2(t)')
% % title('�Ĵ��,Ťת�Ǻͼ��ι���AOA��ʱ��ı仯����')   % �Ĵ��,Ťת�Ǻͼ��ι���AOA��ʱ��ı仯����
% % grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��2�ַ�����⡪������ϵ�£�����Ťת��(x��)����Ťת�������ص��ܺ�M_x
% figure(27)  % ƽ��������������, ת���������ء�ת�������������غ�����������-��Ťת����������
% plot(t/T,M_xtrans,'r-',t/T,M_xrd,'g-',t/T,M_xrotcirc,'m-',t/T,M_xam,'b-',t/T,M_hinge,'c-',t/T,M_inert_x,'y',t/T,M_xtotal2,'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{total,x}  ( including M_{trans,x}  & M_{rd,x}  & M_{rotcirc,x} & M_{am,x} & M_{hinge} & M_{inert,x})  (uN.mm)')
% legend('M_{trans,x}','M_{rd,x}','M_{rotcirc,x}','M_{am,x}','M_{hinge}','M_{inert,x}','M_{total2,x}')
% title('չ��(x-axis)-ƽ��������������,ת����������,ת��������������,����������,�����ظ����غͳ��������������ؼ�������-��ʱ��ı仯����')   % չ������������ʱ��ı仯����
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
% % hold on
% % plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��Ťת��Ĺ���
% % P_xtrans=-M_xtrans.*omega_x*10^-6;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% % P_xrd=-M_xrd.*omega_x*10^-6;  % mW
% % P_xrotcirc=-M_xrotcirc.*omega_x*10^-6;  % mW
% % P_xam=-M_xam.*omega_x*10^-6;  % mW
% % P_inert_x=-M_inert_x.*omega_x*10^-6;  % mW
% % P_totalx=-M_totalx.*omega_x*10^-6;  % mW
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ƽ���������ع���  % uW
% % P_xtrans=-0.1*M_xtrans.*omega_x*10^-3; 
% P_xtrans=-M_xtrans.*omega_x*10^-3;  % uN.mm*rad*s^-1=10^-9N*m*s^-1=10^-9W=10^-6mW=10^-3uW
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % ת���������ع���  % uW
% % P_xrd=-0.5*M_xrd.*omega_x*10^-3;                           % OK
% % M_xrd=-C_RD*omega_x.*abs(omega_x)*M_xrdcoeff*10^(-3);   % ת����������
% P_xrd=-M_xrd.*omega_x*10^-3;                
% % P_xrd=-M_xrd.*omega_x*10^-3; 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % ת���������ع���  % uW  
% % M_xrotcirc(i,1)=k_xRot*C_R*sign(omega_x).*omega_x(i,1).*abs(omega_h(i,1)).*Y_rcpnd_rot(i,1)*M_xRotcoeff*10^(-3); %  ת����������
% % P_xrotcirc=-0.5*(-sign(alpha2(i,1))).*M_xrotcirc.*omega_x*10^-3;  
% % P_xrotcirc=-0.5*M_xrotcirc.*abs(omega_x)*10^-3;    % OK
% P_xrotcirc=-M_xrotcirc.*abs(omega_x)*10^-3; 
% % P_xrotcirc=-M_xrotcirc.*omega_x*10^-3;          % ���������ʽ����ǰ����M_xrotcirc(i,1)��Ҫ����sign(omega_x)*
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % ���������ع��ʺͳ������������ع���  % uW
% P_xam=-M_xam.*omega_x*10^-3;
% P_inert_x=-M_inert_x.*omega_x*10^-3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % M_totalx=M_xtrans+M_xrd+M_xrotcirc+M_xam+M_inert_x;    % ��ƽ��: չ��Ťת��(x��)
% % P_totalx=-M_totalx.*abs(omega_x)*10^-3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_aerox=P_xtrans+P_xrd+P_xrotcirc+P_xam;
% P_totalx=P_xtrans+P_xrd+P_xrotcirc+P_xam+P_inert_x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(28)  % ƽ��������������, ת���������ء�ת�������������غ�����������-��Ťת��Ĺ���
% plot(t/T,P_xtrans,'r-.',t/T,P_xrd,'g--',t/T,P_xrotcirc,'m:',t/T,P_xam,'b-','LineWidth',2)
% hold on
% plot(t/T,P_aerox,'k-','LineWidth',3)
% hold on
% plot([0.9,4.05],[0,0],'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itP_{aero,x}  ( including P_{trans,x}  & P_{rd,x}  & P_{rotcirc,x} & P_{am,x})  (mW)')
% legend('P_{trans,x}','P_{rd,x}','P_{rotcirc,x}','P_{am,x}','P_{aero,x}')
% title('չ��(x-axis)-ƽ��������������,ת����������,ת��������������,���������ؼ�������������Ťת��Ĺ���-��ʱ��ı仯����')   % չ������������ʱ��ı仯����
% grid on
% axis([1,3,-10,10])
% set(gca,'XTick',(1:0.1:3))
% % hold on
% % plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(29)  % �������غ�����������,���������������ؼ���������Ťת��Ĺ���
% plot(t/T,P_aerox,'b-',t/T,P_inert_x,'c:','LineWidth',2)
% hold on
% plot(t/T,P_totalx,'k-','LineWidth',3)
% hold on
% plot([0.9,4.05],[0,0],'k-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itP_{total,x}  ( including P_{aero,x} & P_{inert,x})  (mW)')
% legend('P_{aero,x}','P_{inert,x}','P_{total,x}')
% title('չ��(x-axis)-����������,���������������ؼ���������Ťת��Ĺ���-��ʱ��ı仯����')   % չ������������ʱ��ı仯����
% grid on
% axis([1,3,-10,6.5])
% set(gca,'XTick',(1:0.1:3))
% % hold on
% % plot(t/T,phi*180/pi,'r:',t/T,psi*180/pi,'b:',t/T,alpha2*180/pi,'g:','LineWidth',1.5)  %ת��Ϊms �� ����degree   *10^3   *180/pi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ�ģ�顪������ת���ᡪ���������ط���%%%%%%
%%%%���߲���, �ֱ��ǣ�%%%%%%%%%%%%%%��λ: (mN.mm)��(uN.m)
%%%%��һ���֡���ƽ��������������������%%%%%%%
%%%%�ڶ����֡���ת��������������%%%%%%%%%%%
%%%%�������֡���ת������������%%%%%%%%%%%%
%%%%���Ĳ��֡���ת����������%%%%%%%%%%%%%
%%%%���岿�֡����Ĵ���ת�������ĵ��Իظ�����%%%%
%%%%�������֡����Ĵ�����������%%%%%%%%%%%%
%%%%���߲��֡�����������������%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���� ƽ�������������������������ء����Ƴ�ƽ���µ�������
% ���������������Ϊ��10�� %% ��λ: (mN.mm)��(uN.m)
% Omega1; Omega2; Omega3; Omega4;          �����������������������Գ���ò������
%%%%%%%%%%%%%%%%%%%%%%%%%
% % M_n=(-Omega1*sign(dphi).*dphi.^2.*C_N)*(10^-12*10^6);  % �������ų�ƽ��������(y��)��������bending Motion along chord
% % mg.mm^2*rad*s^-2=(10^-12)kg*m*s^-2*m=10^(-12)N*m=10^-3 uN*mm   
% M_ztrans1=sign(dphi)*Omega1.*dphi.^2.*C_N*(10^-3);  % �������ų�ƽ��������(y��)���������׸���������������ͬ
% k_ztrans=1;  
M_ztrans=k_ztrans*sign(alpha2).*(I1z.*abs(C_N).*omega_h.^2)*10^(-3);   % I1z=0.0162; % ��λ�� mg.mm^2=10^-3uN*mm=10^-6 mN*mm
% % figure(30)
% % plot(t/T,M_ztrans1,'r-',t/T,M_ztrans,'g-','LineWidth',2)  
% % grid on
% % axis([0.9,4.05,-inf,inf])
% r_xcop_tr=M_ztrans./F_ytran;    % ����չ������
% % r_xcop_traver=mean(r_xcop_tr)                 %  r_xcop_traver=-2.3692;
% R_wingeff=3.004;    %��Ч��򳤶�(mm)  
% r_xcopnd_tr=r_xcop_tr/R_wingeff;   
% r_xcopnd_traver=mean(r_xcop_tr)/R_wingeff    % r_xcopnd_traver =-0.7887;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ����� ת�������������������������ء����Ƴ�ƽ���µ�������
% C_R=1.55;    % ��ϵ��������
% M_zrot=I2z*C_R*omega_x.*V_nonr*10^(-3); % V_nonr=omega_h=dphi
C_R=1.55;    % ��ϵ��������
% k_zrot=1;
M_zrot=k_zrot*I2z*C_R*omega_x.*abs(omega_h)*10^(-3);  % I2z=0.0069; ��λ��: mg.mm^2*(rad*s^-1)^2=10^-6mN.mm=10^-3uN.mm
% r_xcop_rot=M_zrot./F_yrot;                          % ����չ������
% r_xcopnd_rot=-abs(r_xcop_rot)/R_wingeff;
% % r_xcop_rotaver=-mean(abs(r_xcop_rot))   % r_xcop_rotaver=-2.1408
% r_xcopnd_rotaver=-mean(abs(r_xcop_rot))/R_wingeff    % r_xcopnd_rotaver=-0.7126;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ƽ��������ת���������������õ�չ��ƽ��ѹ��λ�þ��Ĵ���ľ���
% r_xcopnd=(M_ztrans+M_zrot)./(F_ytran+F_yrot); 
% % r_xcopndaver=mean(r_xcopnd)    %r_xcopndaver=-2.9352(XXX);
% % r_xcopndaver1=(r_xcopnd_traver+r_xcopnd_rotaver)/2;   % r_xcopndaver1=-2.2550;
%% ƽ��������ת����������������չ��ѹ����ʱ��仯����
% figure(31) 
% plot(t/T,r_xcopnd_tr,'b-.',t/T,r_xcopnd_rot,'g-.',t/T,r_xcopnd,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('ƽ��������ת����������������չ��ѹ����Ťת��֮��ľ���  (mm)')
% legend('r_{xcopnd,tr}(t)','r_{xcopnd,rot}(t)','r_{xcopnd}(t)')
% title('ƽ��������ת����������������չ��ѹ����ʱ��ı仯����')
% grid on
% hold on
% % plot(t/T,F_ytran,'r-.',t/T,F_yrot,'g-.',t/T,F_ytran+F_yrot,'k-',t/T,M_ztrans+M_zrot,'k:','LineWidth',2)
% plot(t/T,F_ytran+F_yrot,'k-',t/T,M_ztrans+M_zrot,'k:','LineWidth',2)
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%% �����湥�Ǳ仯��չ��ƽ��������ת��������������ѹ��
% z_cop=[t,alpha2*180/pi,r_xcopnd_tr,r_xcopnd_rot,r_xcopnd]; 
% xlswrite('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\aerodynamic_moment\aerodynamic_moment_model2_2\z_cop.xlsx',z_cop,'sheet1','A1:E2000');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������� �������������������������ء����Ƴ�ƽ���µ������� % I5y��I6y��λ�� mg.mm=*10^(-9) kg.m  %��λ��10^6uN
% % F_yadd1=-(I5y*(domega_z+omega_x.*omega_y)-I6y*domega_x)*10^(-3); % �޸���ԭ���Ƴ��Ĺ�ʽ��������,�Ҳ���/4
% F_yadd1=(-I5y*(domega_z+omega_x.*omega_y)+I6y*domega_x)*10^(-3);
% M_zadd=-(I7z*(domega_z+omega_x.*omega_y)-I_xzam*domega_x)*10^(-3); % �޸���ԭ���Ƴ��Ĺ�ʽ��������,�Ҳ���/4
% k_za=0.35; 
% k_za=1;  % ��ϵ��������ô��
M_zadd=k_za*(-I7z*(domega_z+omega_x.*omega_y)-I_xzam*domega_x)*10^(-3); % ԭ���Ƴ��Ĺ�ʽ������Ϊ����ЩŶ   % ��ʼ��ʱ��(-)
% r_xcopnd_add=M_zadd./F_yadd1;
% % r_xcopnd_addaver=trapz(t,r_xcopnd_add)/(3*T);    % r_xcopnd_addaver=0.7934;
% r_xcopnd_addaver=mean(r_xcopnd_add);   % r_xcopnd_addaver=0.7320;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���Ĳ��� ת���������ء����Ƴ�ƽ���µ�������
% C_RD=5.0;     %ת������������ϵ��:C_RD=C_rd��(3,6); ����C_rd=C_dmax=3.4;
% uN.mm   % (rad/s)^2*mg.mm^2=mg.mm/s^2.mm=10^(-3) uN.mm
% C_RD2=0.05*C_RD;   
% C_RD2=1*C_RD;   % ��ϵ��������ô��
M_zrd=-C_zRD*omega_x.*abs(omega_x)*M_zrdcoeff*10^(-3);   % ��������һ��ʼ��˳ʱ��� % ע���������ȡ����omega_x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ���岿�� �Ƴ�ƽ���µ�����������֮�͡�������������������
% M_ztotal=M_ztrans+M_zrot+M_zadd+M_zrd;
% % M_inert_z=-x_com*F_inert_y*10^3; % ���򡪡���������������;  ��ʼ��ʱ��(-)
M_ztotal=M_ztrans+M_zrot+M_zadd+M_zrd+M_inert_z;  % ��ƽ��: ����(z��)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(32)
% hold on
% plot(t/T,M_ztrans,'r-',t/T,M_zrot,'g-',t/T,M_zadd,'b-',t/T,M_zrd,'c-',...
%        t/T,M_inert_z,'m-',t/T,M_ztotal,'g:','LineWidth',2)  
% xlabel('\itNormalized time')
% ylabel('\itM_{z,trans} & M_{z,rot} & M_{z,add} & M_{z,rd} & M_{inert,z}(t) & M_{total,z}  (uN.mm)')
% legend('M_{z,trans}(t)','M_{z,rot}(t)','M_{z,add}(t)','M_{z,rd}(t)','M_{inert,z}(t)','M_{total,z}(t)')
% title('����(������ϵz��)����������ʱ��ı仯����')   % ��������������ʱ��ı仯����
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ģ�顪������ת���ᡪ���������ط���%%%%%
%%%%��������, �ֱ��ǣ�%%%%%%%%%%%%%%��λ: (mN.mm)��(uN.m)
%%%%��һ���֡���ƽ��������������������%%%%%%
%%%%�ڶ����֡���ת��������������������%%%%%%
%%%%�������֡�����������������%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���֡���ƽ�������������������ء���������
% ���������������Ϊ��10�� %% ��λ: (mN.mm)��(uN.m)
% Omega1; Omega2; Omega3; Omega4;         �����������������������Գ���ò������
%%%%%%%%%%%%%%%%%%%%%%%%%
% M_t=(-Omega1*sign(dphi).*dphi.^2.*C_T)*(10^-12*10^6);    % �������ų�ƽ�淨����(z��)��������Plunging Motion
% mg.mm^2*rad*s^-2=(10^-12)kg*m*s^-2*m=10^-12N*m=10^-3 uN*mm
% M_tang=-sign(dphi)*Omega1.*dphi.^2.*C_T*(10^-3);    % uN*mm %  �������ų�ƽ�淨����(z��)���������׸���
% M_ytrans1=M_tang;            % ��ƽ��: ����(y��)
M_ytrans=-C_Ty1*sign(alpha2).*(I1z.*C_T.*omega_h.^2)*10^(-3);   % I1z=0.0162; % ��λ�� mg.mm^2
% figure(33)
% plot(t/T,M_ytrans1,'r-',t/T,M_ytrans,'g-','LineWidth',2)  
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ����֡���ת�������������������ء���������
% C_R=1.55;    % ��ϵ�������⡪��������������XXXXXXXXXXXXXX
% M_zrot=I2z*C_R*omega_x.*V_nonr*10^(-3);   % V_nonr=omega_h=dphi
% I2z=0.0069; ��λ��: mg.mm^2  * (rad*s^-1)^2=mN.mm=10^(-3)uN.mm
M_yrot=C_Ty2*I2z*C_T.*omega_x.*abs(omega_h)*10^(-3);    % C_Tȫ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������֡����������������ء���������
% M_ytotal=M_ytrans+M_yrot;
M_ytotal=M_ytrans+M_yrot+M_inert_y;  % M_inert_y=x_com*F_inert_z*10^3;  % ���򡪡���������������;  ��ʼ˳ʱ��(-)
% figure(34)
% plot(t/T,M_ytrans,'r-',t/T,M_yrot,'b-',t/T,M_inert_y,'g-',t/T,M_ytotal,'k-','LineWidth',2)  
% xlabel('\itNormalized time')
% ylabel('\itM_{trans,y} & M_{rot,y} & M_{inert,y} & M_{total,y} (uN.mm)')
% legend('M_{trans,y}(t)','M_{rot,y}(t)','M_{inert,y}(t)','M_{total,y}(t)')
% title('����(������ϵy��)����������ʱ��ı仯����')   % ��������������ʱ��ı仯����
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ģ�顪����������ϵ�µĸ��᷽���������ط���%%%%%
%%%%��������, �ֱ��ǣ�%%%%%%%%%%%%%%��λ: (mN.mm)��(uN.m)
%%%%��һ���֡�����������ϵ�µ�X�᷽������%%%%%%
%%%%�ڶ����֡�����������ϵ�µ�Y�᷽������%%%%%%
%%%%�������֡�����������ϵ�µ�Z�᷽������%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81;         % �����������ٶ�:g=9.821N/kg=9.821*10^6/10^6=9.821 uN/mg  ����g�Ĺ��ʵ�λ��m*s^-2��N*kg^-1
M_body =1.8;  %��Ӭ���������(mg)��Science����M_body =1.8e-06;(kg)
W=M_body*g;
R_wingeff=3.004;    %��Ч��򳤶�(mm)��
norm_denom=W*R_wingeff;    % ��λ:  uN.mm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ��һ���֡�����������ϵ�µ�X�᷽������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_X=cos(phi).*M_xtotal+sin(phi).*cos(psi).*M_ytotal-sin(phi).*sin(psi).*M_ztotal;
% figure(35)
% plot(t/T,2*M_X/norm_denom,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{X} (uN.mm)')
% legend('M_{X,pitch}')
% title('��������ϵ�µ�X�᷽��pitch������-�������')  % ������������
% grid on
% % axis([0.9,4.05,-inf,inf])
% axis([0.9,2.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ����֡�����������ϵ�µ�Y�᷽������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_Y=-sin(phi).*M_xtotal+cos(phi).*cos(psi).*M_ytotal-cos(phi).*sin(psi).*M_ztotal;
% figure(36)
% plot(t/T,-M_Y/norm_denom,'r-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{Y} (uN.mm)')
% legend('M_{Y,roll}')
% title('��������ϵ�µ�Y�᷽�����ء�roll������-�������')   % ������������
% grid on
% % axis([0.9,4.05,-inf,inf])
% axis([0.9,2.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �������֡�����������ϵ�µ�Z�᷽������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1)�����Ĵ���ת�����������򡪽����ظ�����
% T_transRat=3000*10^-3;          % ������    % ��λ: rad/mm; 
T_transRat=2857*10^-3;              % ������    % ��λ: rad/mm; 
k_act=300;                                    % ��λ:  N/m=mN/mm  % 2011-ICIRS-System identification
% �ο�: 2010-BB-Distributed power and control actuation
% k_zh=k_act/T_transRat^2;          % k_zh =33.3333=mN*mm^-1* rad^-2*mm^2=mN.mm* rad^-2 
k_eq=344.8;                                   % ��λ:  N/m=mN/mm
% k_zh1=2*10^(-3);
k_zh1=0.15*10^(-3);
k_zh=k_zh1*(k_eq-k_act)/T_transRat^2;   % k_trans =5.4885; % mN.mm/rad;
% k_zh=5.4872;                              % k_trans=5.4872;  % uN.m/rad=10^-3*10^3 mN.mm/rad;
M_zhinge=k_zh*phi*10^(3);         % ��λ��: mN.mm=*10^(3)uN.mm
% figure(37)
% plot(t/T,M_zhinge,'c-')  
% xlabel('\itNormalized time')
% ylabel('\itM_{z,hinge}   (uN.mm)')      
% legend('M_{z,hinge}(t)')
% title('�����ظ�����(չ��:x-axis)��ʱ��ı仯����')
% grid on
% axis([0.9,4.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2)�����Ĵ����������������
% F_act=F_0*sin(2*pi*f*t );  % �ο�: 2010-BB-Distributed power and control actuation
T_transRat=3000*10^-3;    % ������    % ��λ: rad/mm;        % T_transRat=3000rad/um;  % 2012-BB-��ʽ��-Conceptual design 
% F_act=55.9;                         %��λ: mN % for phi_halfmax =57.5823��;                           % F_b = 0.0559N   0.0594N 
% T_transRat=3300*10^-3;     % ������    % ��λ: rad/mm         % 2008-IEEE TR-RJ Wood
% F_act=136/2;                        % F_act=123/2; % ��λ: mN
% F_act=112;      % ��λ: mN   % ��Ҫ��� PHI_pp=1.1487*180/pi=65.82��;
F_act=56;          % ��λ: mN  % ƥ���������������0.3mN.mm
M_ampl=F_act/T_transRat;
w =1185.6;
T=2*pi/w;    % f=1/T;   % f=188.7;
% k_zact=0.25*10^(2);
% k_zact=6.75*10^(3);
k_zact=5*10^(-3);
% k_zact=7.5*10^(2);
% delta_act=0; 
delta_act=-1.571;
% t_range=linspace(0.0052824335,0.0052824335+5*T,1000);
M_zact=k_zact*M_ampl*cos(w*t+delta_act)*10^(3);  % ��λ��: mN.mm=*10^(3)uN.mm
% figure(38)
% plot(t/T,M_zact,'r-')
% xlabel('\itNormalized time');  ylabel('\itM_{z,act} (mN.mm)');
% legend('M_{z,act}(t)')
% title('�Ĵ����������������')
% grid on
% axis([0.9,4.05,-inf,inf])
% set(gca,'XTick',(0.9:0.1:4.05))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_Z=sin(psi).*M_ytotal+cos(psi).*M_ztotal;
M_Z_ActHinge=M_Z-M_zact+M_zhinge; 
% figure(39)
% plot(t/T,M_Z/norm_denom,'r-',t/T,M_Z_ActHinge/norm_denom,'b-','LineWidth',2)
% xlabel('\itNormalized time')
% ylabel('\itM_{Z} & M_{Z,ActHinge} (uN.mm)')
% legend('M_{Z,yaw}','M_{Z,ActHinge,yaw}')
% title('��������ϵ�µ�Z�᷽��yaw������-�������')  % ������������
% grid on
% % axis([0.9,4.05,-inf,inf])
% axis([0.9,2.05,-inf,inf])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����science��е��Ӭ����������������ݡ�����һ��֮�������
%%����Ķ�����������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % output_3=[t_NOfreq,My_norm_stroke_butterfilt_mean,My_norm_rotation_butterfilt_mean,My_norm_deviation_butterfilt_mean];
% pitchmoment_science=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\aerodynamic_moment\aerodynamic_moment_model2_2\force_pitch_y_moment_ForceModul.xlsx','A1:D1145'); % ��������
%%�������������ȷ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output_2=[t_NOfreq,Mx_norm_all_butterfilt_steady,My_norm_all_butterfilt_steady,Mz_norm_all_butterfilt_steady];
pitchmoment_science=xlsread('D:\KXJ\PassiveRot_dynamic_Science_fruitfly\wing_parameter\datanalysis_science_fruitfly\robotForcesTorques\ForceModulations\aerotorque_for_steady_wingbeat.xlsx','A1:D1145');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_NOfreq1=pitchmoment_science(:,1);  % t_NOfreq(1,1)=0;
% Mx_norm1=pitchmoment_science(:,2); % ��������
My_norm1=pitchmoment_science(:,3);     % ��������
% Mz_norm1=pitchmoment_science(:,4); % ƫ������
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_NOfreq=[t_NOfreq1+0.0052824335;t_NOfreq1+T+0.0052824335;t_NOfreq1+2*T+0.0052824335];
My_norm=[My_norm1;My_norm1;My_norm1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%��������ϵ�µ��������ط�����ʱ��ı仯������ʵ����Խ���ĶԱ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(43)
plot(t_NOfreq/T,2*M_X/norm_denom,'r--',t_NOfreq/T,My_norm-0.02,'b-','LineWidth',2)
% % 1.525�������ǳ�������  % 1.65���������ǳ������������ϸ������Ͻ�Ӧ��������������������������
xlabel('\itNormalized time')
ylabel('\itM_{X,norm,cal} & M_{y,norm,str} (uN)')
legend('M_{X,norm,cal}','M_{y,norm,str,exp}')
title('����������,׼��̬����Ԥ��ĸ������غ�ʵ���õĳ�̽�����ĸ������صĶԱ�')   
grid on
axis([0.9,4.05,-inf,inf])
set(gca,'XTick',(0.9:0.1:4.05))
hold on
L=length(t);
plot([0,t(L)/T],[0,0],'k-','LineWidth',2);     %��x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ĸ���������ϵ���ĳͷ����Լ��
coeff_con=aeroM_coeff_constraint(x);
% ������������ϵ���ĳͷ����penaltyfun
s=2000;
penaltyfun=s*coeff_con;           % penaltyfun =;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdev_pitch=std(((2*M_X/norm_denom)-(My_norm-0.02)),0,1);   % ����S1(0)������(1)Ԫ�صı�׼����
obj_function=stdev_pitch+penaltyfun;   % obj_function=;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc  % Elapsed time is 218.831766 seconds.
