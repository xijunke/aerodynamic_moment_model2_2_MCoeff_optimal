function Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha)
% �޸�ʱ�䡪��2014��12��21��,23:36
%ת����ƫ�������е��ƫ��������,��Ťת������ƫ��C_maxy֮��
%% (1) �����ת���������ء�����⾻ѹ�ĵ�������λ��Y_cpnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R=16.0148-0.88=15.1348;
% �� r_nd��(0.88/15.1348,13.4427/15.1348)
% clear all; clc;
% C_avereff=0.8854;          % mm    
R_wingeff=3.004;          %��Ч��򳤶�(mm) 
xr=0.3289;                     % x-root offset  \mm
xr_nd=xr/R_wingeff;      % x-root offset  ������չ��ƫ�þ���
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ���-ƫ��-����ϵԭ��ľ���
% R_proximal=xr;                                                    % xr=3.19;     %RJ Wood��Ƶĳ��\mm
% R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood��Ƶĳ��\mm
% x=linspace(R_proximal,R_distal,200);
% x_mod_Root=0.636;                                            % mm
% C_maxy=0.138924474377504;  % ��Գ���wing_model_88_yaxis��: ��122��; C_maxy =0.1389; 
% yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root-C_maxy;  
% yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root-C_maxy;
% % ����(1)������ǰ��Ե�����������ٻ�yr_leadnd����yr_trailnd�������
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd-0.156013;
% Cr_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms r_nd
% �����ǳ�ǰԵ�����������Ťת���λ�ò�ͬ��Ҫ���зֶκ�������
yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% �������ҳ��ֲ�Ϊ6�׶���ʽ
C_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();        %���ú���kenimatics_wing_and_AoA
% % size(wing_kenimatics)     %  (1000*12)
% alpha=wing_kenimatics(:,4);
% alpha=pi/2.4;                                  % ��λ��rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ����ѹ�ĵ�ķֲ�������Ťת�����۵ķ���ı�
d_cprnd=0.82*abs(alpha)/pi+0.05;     % ��ת���������ص�����ѹ��λ��; %��abs(alpha)��(pi/4,pi/2)ʱ��d_cprnd��(0.255,0.46);
% yr_cpnd=yr_nd+yr_leadnd-C_nd*d_cprnd; 
%% ��һ����������������Ŷ
%  if d_cprnd>0.25  % ��������d_cprnd������,������һ�������������(yr_leadnd-yr_nd)/C_nd=0.25;
%      yr_cpnd=C_nd*d_cprnd-yr_leadnd;           %  z3=-(d_cp-C_025);
% % ���������⣬�������ĳһ���㶨����ʱ,����չ��ѹ�ĵķֲ�������Ťת��ǰ����߿����ں��棬����Ӧ�÷�������л��ֵ���
%  elseif d_cprnd<=0.25;   % ��������d_cprnd��������������һ�������������(yr_leadnd-yr_nd)/C_nd=0.25;
%      yr_cpnd=yr_leadnd-C_nd*d_cprnd;           %  z3=C_025-d_cp;
% % ���������⣬�������ĳһ���㶨����ʱ������չ��ѹ�ĵķֲ�������Ťת��ǰ����߿����ں��棬����Ӧ�÷�������л��ֵ���
%  end
%% �ڶ������������Կ���ѡȡ��Ӧ��չ������ѹ��λ������Ӧ��Ƭ�����з�����������ѹ�� 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) ѡȡ��Ӧ��չ������ѹ��λ������Ӧ��Ƭ�����з�����������ѹ��                                  
% % r_xcopnd_rot =0.7126;                               % r_xcopnd_rot=I2/F_ndRot;  % *R_wingeff 
% % R_rot=(R_wingeff+xr)*r_xcopnd_rot;
% % r_nd_rot=(R_rot-xr)/R_wingeff;                  % r_nd_rot =0.6811;
r_nd_rot=0.712638285096407;                          % �� 0.7126; 
% Ratio_axis=vpa(yr_leadnd./C_nd, 5);               % vpa(f,d); d��ָ��Ч���ֵ�λ��
% Ratio=inline(vectorize(Ratio_axis),'r_nd');
% Ratio_axisrot=Ratio(r_nd_rot)                         % Ratio_axisrot =0.2875; ��������һ���ǳ���Ҫ�Ĳ���
%% ȷ��ƽ��������չ��ѹ�Ĵ���ӦƬ����������ǰԵ���������ҳ� @ r_nd_rot =0.6811;
yr_leadnd_rot1=inline(vectorize(yr_leadnd),'r_nd');
yr_leadnd_rot=yr_leadnd_rot1(r_nd_rot);         % yr_leadnd_rot =0.5221
C_nd_rot1=inline(vectorize(C_nd),'r_nd');
C_nd_rot=C_nd_rot1(r_nd_rot);    % C_nd_rot=1.2877=C_cop_rot =1.1401/0.8854=1.28767;  
% alpha_crit_abs=(Ratio_axis-0.05)*pi/0.82*180/pi  % ��d_cprnd=Ratio_axisʱ, ȷ���ٽ����ֵ�Ƕ�alpha_crit_abs=81.2374;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if d_cprnd>Ratio_axis  % ��������d_cprnd��������������һ����������0.25
%      yr_cpnd=-(C_nd_rot*d_cprnd-yr_leadnd_rot);     %  z3=-(d_cp-C_025);   % ����  %��abs(alpha)��(pi/4,pi/2)ʱ��d_cprnd��(0.255,0.46);
% % ���������⣬�������ĳһ���㶨����ʱ������չ��ѹ�ĵķֲ�������Ťת��ǰ����߿����ں��棬����Ӧ�÷�������л��ֵ���
% elseif d_cprnd<=Ratio_axis
%      yr_cpnd=yr_leadnd_rot-C_nd_rot*d_cprnd;    %  z3=C_025-d_cp;            % ����  %��abs(alpha)��(pi/4,pi/2)ʱ��d_cprnd��(0.255,0.46);
% % ���������⣬�������ĳһ���㶨����ʱ������չ��ѹ�ĵķֲ�������Ťת��ǰ����߿����ں��棬����Ӧ�÷�������л��ֵ���
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr_cpnd=yr_leadnd_rot-C_nd_rot*d_cprnd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx6=(r_nd+xr_nd)*C_nd.^2;                  % ������������F_ndRot��ԭʼ��������
fx7=expand(fx6);
F_ndRot=double(int(fx7,r_nd,0,1));       % Result: F_ndRot=0.74851
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ����������������(nondimention_aerodynamic_component)����⡪��F_ndRot
% ע�⡪���öγ����мǲ����޸ģ�ǰ��ֻҪ��֤������ȷ���������ҳ��ֲ����ɡ�
%���µĹ�ʽӦʹ�ú����������ٵ��ҳ��ֲ���ʽC_nd
% R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1));   % ��������صĻ�ת�뾶��ƽ��
% R1nd1=double(int(r_nd*C_nd,r_nd,0,1));        % һ������صĻ�ת�뾶
% % S_nd=double(int(C_nd,r_nd,0,1));               % �����ٳ����
% F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;      % ʹ����������Ҳ��ȷ; ���:F_nd2 =0.5024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd_rot=int(yr_cpnd*fx7,r_nd,0,1)/F_ndRot;
% disp(['��ѹ�ĵ�������λ��Y_cpnd_rot(alpha)=' num2str(Y_rcpnd_rot)  ' ���ٵ�λ������mm'])
% Y_rcpnd_rot=abs(double(int(yr_cpnd*fx7,r_nd,0,1))/F_ndRot);
Y_rcpnd_rot=double(int(yr_cpnd*fx7,r_nd,0,1))/F_ndRot;
