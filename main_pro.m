clc
clear
%�������������

global t %��¯ʱ��
global de %��������1�ľ���
global de1 %��������2�ľ���
%  d1=1;d2=1;

v=7/6;%�ٶ�
fujian=xlsread('����.xlsx');
plot(fujian(:,1),fujian(:,2));%���Ƹ�������ͼ��

%���������Է���ʱȡ��ע��
% for de=1.5:0.1:2.5
% for de1=60:70
F=zeros(size(fujian,1),1);
delta=0.5/0.02; %ÿ��0.5s���Ƶļ��
k1=1;k2=1;k3=1;

%��a��k1/k���д��·�Χ�Ĳ����������
for aaa=0.000015:0.000001:0.000025
    for kkk=3000:100:8000
        zzz=chuanrefangcheng(kkk,aaa);
        t=t-t(1,1);
        n=find(t==19);
        %         for i=1:size(zzz,1)
        %            if(zzz(i,16)>30+273)
        %               n=i;
        %                break;
        %            end
        %         end
        for i=n:delta:size(zzz,1)
            xx(k1,k2)=zzz(i,16);
            FF(k1,k2)=(fujian(k1,2)+273-xx(k1,k2)).^2;%��С����
            k1=k1+1;
        end
        fengzhi(k2)=max(xx(:,k2));
        k1=1;
        k2=k2+1;
    end
    FF=sum(FF);
    xxx{k3}=xx;
    FFF{1,k3}=FF;
    FFF{2,k3}=fengzhi;
    k3=k3+1;
    k2=1;
end
for i=1:size(FFF,2)
    for j=1:size(FFF{1,1},2)
        if(abs(FFF{2,i}(j)-max(fujian(:,2)+273))>2)
            FFF{1,i}(j)=1000000;
        end
    end
end

%ѡ����С�������Ž�
for i=1:size(xxx,2)
    [fm(i),ind(i)]=min(FFF{1,i});
    [fmin,indx]=min(fm);
end

%���������Է���ʱȡ��ע��
%         error(d1,d2)=fmin;
%         d1=d1+1;
%     end
%     d1=1;
%     d2=d2+1;
% end

%������Ͻ��
figure;
plot(fujian(:,1),xxx{indx}(:,ind(indx))-273,'Linewidth',1.2);
hold on;
plot(fujian(:,1),fujian(:,2),'Linewidth',1.2);
grid on;
ylabel('���Ӱ����������¶�/��');
xlabel('ʱ��/s');
title('¯��������С�����������ͼ');
legend('��Ͻ��','��������');

%���������Ͻ��
aaa=0.000015:0.000001:0.000025;
kkk=3000:100:8000;

disp('������С���˷���ϵõ��Ĳ���aΪ��');
disp(aaa(indx));
disp('������С���˷���ϵõ��Ĳ���k1/kΪ��')
disp(kkk(ind(indx)));

%�����Է���ͼ��
% mesh(de,de1,error,'Linewidth',1.2);
% title('�¶ȹ������������ȷ������');
% xlabel('��ȴ��֮ǰ���¶ȹ����������/m','Rotation',0);
% ylabel('��ȴ��֮����¶ȹ����������/m');
% zlabel('���ƽ����');