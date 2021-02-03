clc
clear
%��������� �������

%��ֵ�����½� ���г�ֵʹ�õ����ʵ������
x0=[180.59,196.14,226.68,264.32,85.9];
lb=[165,185,225,245,65];
ub=[185,205,245,265,100];
ms = MultiStart('TolX',1.0e-10,'MaxTime',300) ;%
problem = createOptimProblem('fmincon', 'objective',@func_pro4,'x0',x0, 'lb',lb,'ub',ub);
[a,fminm,flagm,outptm,manyminsm] = run(ms,problem,60);%���

%����鿴���
for i=1:length(manyminsm)
    result(i,1)=manyminsm(i).Fval;
    result(i,2)=manyminsm(i).Exitflag;
    Tvv{i}=manyminsm(i).X;
end

