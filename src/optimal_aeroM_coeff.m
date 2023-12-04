%% optimal_aeroM_coeff
%% ��һ��������Call ga�����Ŵ��㷨
tic
matlabpool open 8
% function obj_function=objfunction_aeroM_coeff(x)
ObjFunction=@objfunction_aeroM_stdev;  % fitness and constraint functions��������Ŀ�꺯��objfunction_aeroM_coeff
nvars=14;   % Number of variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=gaoptimset('Display','diagnose','Display','iter','StallTimeLimit',7200,'PopulationSize',100,...
               'PlotFcns',{@gaplotbestf,@gaplotexpectation,@gaplotstopping,@gaplotbestindiv,@gaplotscores,@gaplotrange},...
               'UseParallel','always'); %���м��㡪����
% options=gaoptimset('PopulationSize',100,'PlotFcns',{@gaplotbestf,@gaplotmaxconstr},'UseParallel','always'); %���м������ % 20150407
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fminsearchoptions = optimset('Display','iter','MaxFunEvals',1000,'MaxIter',1000,'TolCon',1e-6,'TolFun',1e-6,'UseParallel','always');
% 'Display','notify-detailed'
fminsearchoptions = optimset('Display','iter-detailed','MaxFunEvals',1000,'MaxIter',1000,'TolCon',1e-15,'TolFun',1e-15,'UseParallel','always');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options =gaoptimset(options,'Hybridfcn',{@fminsearch,fminsearchoptions});
%%Call ga�����Ŵ��㷨
[x,fval,exitflag,output,population,scores]=ga(ObjFunction,nvars,[],[],[],[],[],[],[],options)
matlabpool close
toc
%%
% output_1=[population,scores];
% xlswrite('D:\KXJ\hybrid_GA_fminsearch_WingM4_4\population_scores.xlsx',output_1,'sheet1','A1:J100');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% �ڶ��������� Call patternsearch����ģʽ����
tic
matlabpool open 8
ObjFunction=@objfunction_aeroM_coeff;  % fitness and constraint functions��������Ŀ�꺯��objfunction_aeroM_coeff
x0=[1, 1, 1, 1, 1.11, 1, 1, 1, 1, 1, 1, 1, 1, 1];      % stdev =;
% psoptimset����Create pattern search options structure
% % Complete poll around current iterate; 'off'
options=psoptimset('CompletePoll','on','CompleteSearch','on','Display','iter','UseParallel','always'); 
[x,fval,exitflag,output]=patternsearch(ObjFunction,x0,[],[],[],[],[],[],[],options)
matlabpool close
toc
%%
% output_1=[population,scores];
% xlswrite('D:\KXJ\hybrid_GA_fminsearch_WingM4_4\population_scores.xlsx',output_1,'sheet1','A1:J100');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


