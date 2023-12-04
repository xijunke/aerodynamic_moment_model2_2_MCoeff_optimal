function coeff_con=aeroM_coeff_constraint(x)
% 建立气动力矩系数约束之外的惩罚函数
% k_xaero=x(1);    % [0,2];
% k_xRot=x(2);      % [0,1.5];
% C_xRD=x(3);      % [0,1.5];
% k_am_x=x(4);     % [0,2.35];
% k_hinge_x=x(5); % [0,10];
% k_inert_x=x(6);   % [0,2];
% k_inert_y=x(7);   % [0,2];
% k_inert_z=x(8);   % [0,2];
% k_ztrans=x(9);    % [0,2];
% k_zrot=x(10);     % [0,1.5];
% k_za=x(11);        % [0,2.35];
% C_zRD=x(12);     % [0,1.5];
% C_Ty1=x(13);     % [0,2];
% C_Ty2=x(14);     % [0,2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LB=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];                                            % Lower bound     % 14个变量 
UB=[2, 1.5, 1.5, 2.35, 10,2, 2,2, 2, 1.5, 2.35, 1.5, 2,2];                     % Upper bound     % 14个变量
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(x);
con_min=LB;
con_max=UB;
zeta=zeros(N,1);
% y=zeros(N,1);
parfor i=1:N
    if x(i) < con_min(i)
        zeta(i)=abs(con_min(i)-x(i))/(con_max(i)-con_min(i));
    elseif x(i) > con_max(i)  % 
        zeta(i)=abs(x(i)-con_max(i))/(con_max(i)-con_min(i));
    else % x(i)>=con_min(i) && x(i)<=con_max(i);  % 无需这句表达式
        zeta(i)=0;  % disp('变量未超出边界约束');  
    end
end
coeff_con=sum(zeta);
end

