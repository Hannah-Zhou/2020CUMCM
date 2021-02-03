clc
clear
%问题三求解 需要进行多次求解

%设定初值和上下界
x0=[180,195,225,265,85];
lb=[165,185,225,245,65];
ub=[185,205,245,265,100];

%设置多起点求解和fmincon求解器
ms = MultiStart('TolX',1.0e-10,'MaxTime',300) ;
problem = createOptimProblem('fmincon', 'objective',@func_pro3,'x0',x0, 'lb',lb,'ub',ub);
[a,fminm,flagm,outptm,manyminsm] = run(ms,problem,40);%求解

%方便查看结果
for i=1:length(manyminsm)
    result(i,1)=manyminsm(i).Fval;
    result(i,2)=manyminsm(i).Exitflag;
    Tvv{i}=manyminsm(i).X;
end