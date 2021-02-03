clc
clear
%��������� ��Ҫ���ж�����

%�趨��ֵ�����½�
x0=[180,195,225,265,85];
lb=[165,185,225,245,65];
ub=[185,205,245,265,100];

%���ö��������fmincon�����
ms = MultiStart('TolX',1.0e-10,'MaxTime',300) ;
problem = createOptimProblem('fmincon', 'objective',@func_pro3,'x0',x0, 'lb',lb,'ub',ub);
[a,fminm,flagm,outptm,manyminsm] = run(ms,problem,40);%���

%����鿴���
for i=1:length(manyminsm)
    result(i,1)=manyminsm(i).Fval;
    result(i,2)=manyminsm(i).Exitflag;
    Tvv{i}=manyminsm(i).X;
end