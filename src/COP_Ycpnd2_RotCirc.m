function Y_rcpnd_rot=COP_Ycpnd2_RotCirc(alpha)
% 修改时间——2014年12月21日,23:36
%转动轴偏离弦向中点的偏移量坐标,在扭转轴向上偏移C_maxy之后
%% (1) 针对旋转轴气动力矩——求解净压心的无量纲位置Y_cpnd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R=16.0148-0.88=15.1348;
% 当 r_nd∈(0.88/15.1348,13.4427/15.1348)
% clear all; clc;
% C_avereff=0.8854;          % mm    
R_wingeff=3.004;          %有效翅膀长度(mm) 
xr=0.3289;                     % x-root offset  \mm
xr_nd=xr/R_wingeff;      % x-root offset  无量纲展向偏置距离
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 翅根-偏离-坐标系原点的距离
% R_proximal=xr;                                                    % xr=3.19;     %RJ Wood设计的翅膀—\mm
% R_distal=R_wingeff+xr;                                        % yr=0.73;    %RJ Wood设计的翅膀—\mm
% x=linspace(R_proximal,R_distal,200);
% x_mod_Root=0.636;                                            % mm
% C_maxy=0.138924474377504;  % 针对程序wing_model_88_yaxis有: 第122行; C_maxy =0.1389; 
% yr_lead=-0.08249*x.^6+0.9167*x.^5-4.04*x.^4+8.872*x.^3-10.06*x.^2+5.674*x-0.413-x_mod_Root-C_maxy;  
% yr_trail=-0.0333*x.^6+0.504*x.^5-2.795*x.^4+7.258*x.^3-8.769*x.^2+3.739*x+0.1282-x_mod_Root-C_maxy;
% % 方案(1)——由前后缘函数的无量纲化yr_leadnd——yr_trailnd——求得
% yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% yr_trailnd =-27.6379*r_nd^6+121.093*r_nd^5-185.804*r_nd^4+125.603*r_nd^3-33.1053*r_nd^2-0.145527*r_nd-0.156013;
% Cr_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms r_nd
% 下面是翅前缘函数——针对扭转轴的位置不同需要进行分段函数处理
yr_leadnd =-68.4639*r_nd^6+208.297*r_nd^5-245.231*r_nd^4+137.467*r_nd^3-36.8594*r_nd^2+4.79235*r_nd-0.156071;
% 无量纲弦长分布为6阶多项式
C_nd =-40.826*r_nd^6+87.204*r_nd^5-59.4267*r_nd^4+11.8645*r_nd^3-3.75408*r_nd^2+4.93788*r_nd-0.0000578215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wing_kenimatics=kenimatics_wing_and_AoA_fruitfly_exp();        %调用函数kenimatics_wing_and_AoA
% % size(wing_kenimatics)     %  (1000*12)
% alpha=wing_kenimatics(:,4);
% alpha=pi/2.4;                                  % 单位是rad
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 弦向压心点的分布—决定扭转轴力臂的方向改变
d_cprnd=0.82*abs(alpha)/pi+0.05;     % 旋转轴气动力矩的弦向压心位置; %当abs(alpha)∈(pi/4,pi/2)时，d_cprnd∈(0.255,0.46);
% yr_cpnd=yr_nd+yr_leadnd-C_nd*d_cprnd; 
%% 第一方案——不够合理哦
%  if d_cprnd>0.25  % 问题在于d_cprnd是向量,而不是一个数，例如假设(yr_leadnd-yr_nd)/C_nd=0.25;
%      yr_cpnd=C_nd*d_cprnd-yr_leadnd;           %  z3=-(d_cp-C_025);
% % 这里有问题，即便针对某一个恒定攻角时,沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
%  elseif d_cprnd<=0.25;   % 问题在于d_cprnd是向量，而不是一个数，例如假设(yr_leadnd-yr_nd)/C_nd=0.25;
%      yr_cpnd=yr_leadnd-C_nd*d_cprnd;           %  z3=C_025-d_cp;
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
%  end
%% 第二方案——可以考虑选取相应的展向气动压心位置所对应的片条进行分析——弦向压心 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) 选取相应的展向气动压心位置所对应的片条进行分析——弦向压心                                  
% % r_xcopnd_rot =0.7126;                               % r_xcopnd_rot=I2/F_ndRot;  % *R_wingeff 
% % R_rot=(R_wingeff+xr)*r_xcopnd_rot;
% % r_nd_rot=(R_rot-xr)/R_wingeff;                  % r_nd_rot =0.6811;
r_nd_rot=0.712638285096407;                          % ≈ 0.7126; 
% Ratio_axis=vpa(yr_leadnd./C_nd, 5);               % vpa(f,d); d是指有效数字的位数
% Ratio=inline(vectorize(Ratio_axis),'r_nd');
% Ratio_axisrot=Ratio(r_nd_rot)                         % Ratio_axisrot =0.2875; ————一个非常重要的参数
%% 确定平动气动力展向压心处对应片条的无量纲前缘和无量纲弦长 @ r_nd_rot =0.6811;
yr_leadnd_rot1=inline(vectorize(yr_leadnd),'r_nd');
yr_leadnd_rot=yr_leadnd_rot1(r_nd_rot);         % yr_leadnd_rot =0.5221
C_nd_rot1=inline(vectorize(C_nd),'r_nd');
C_nd_rot=C_nd_rot1(r_nd_rot);    % C_nd_rot=1.2877=C_cop_rot =1.1401/0.8854=1.28767;  
% alpha_crit_abs=(Ratio_axis-0.05)*pi/0.82*180/pi  % 当d_cprnd=Ratio_axis时, 确定临界绝对值角度alpha_crit_abs=81.2374;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if d_cprnd>Ratio_axis  % 问题在于d_cprnd是向量，而不是一个数，例如0.25
%      yr_cpnd=-(C_nd_rot*d_cprnd-yr_leadnd_rot);     %  z3=-(d_cp-C_025);   % 负号  %当abs(alpha)∈(pi/4,pi/2)时，d_cprnd∈(0.255,0.46);
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
% elseif d_cprnd<=Ratio_axis
%      yr_cpnd=yr_leadnd_rot-C_nd_rot*d_cprnd;    %  z3=C_025-d_cp;            % 正号  %当abs(alpha)∈(pi/4,pi/2)时，d_cprnd∈(0.255,0.46);
% % 这里有问题，即便针对某一个恒定攻角时，沿着展向，压心的分布可能在扭转轴前面或者可能在后面，所以应该分区间进行积分叠加
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yr_cpnd=yr_leadnd_rot-C_nd_rot*d_cprnd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fx6=(r_nd+xr_nd)*C_nd.^2;                  % 无量纲气动力F_ndRot的原始被积函数
fx7=expand(fx6);
F_ndRot=double(int(fx7,r_nd,0,1));       % Result: F_ndRot=0.74851
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 无量纲气动力分量(nondimention_aerodynamic_component)的求解——F_ndRot
% 注意——该段程序切记不得修改，前提只要保证输入正确的无量纲弦长分布即可。
%以下的公式应使用合理的无量纲的弦长分布公式C_nd
% R2nd2=double(int(r_nd^2*C_nd,r_nd,0,1));   % 二阶面积矩的回转半径的平方
% R1nd1=double(int(r_nd*C_nd,r_nd,0,1));        % 一阶面积矩的回转半径
% % S_nd=double(int(C_nd,r_nd,0,1));               % 无量纲翅面积
% F_nd2=R2nd2+2*xr_nd*R1nd1+xr_nd^2;      % 使用这句计算结果也正确; 输出:F_nd2 =0.5024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y_rcpnd_rot=int(yr_cpnd*fx7,r_nd,0,1)/F_ndRot;
% disp(['净压心的无量纲位置Y_cpnd_rot(alpha)=' num2str(Y_rcpnd_rot)  ' 量纲单位可能是mm'])
% Y_rcpnd_rot=abs(double(int(yr_cpnd*fx7,r_nd,0,1))/F_ndRot);
Y_rcpnd_rot=double(int(yr_cpnd*fx7,r_nd,0,1))/F_ndRot;

