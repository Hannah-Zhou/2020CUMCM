clc
clear
%参数求解主函数

global t %过炉时间
global de %过渡区域1的距离
global de1 %过渡区域2的距离
%  d1=1;d2=1;

v=7/6;%速度
fujian=xlsread('附件.xlsx');
plot(fujian(:,1),fujian(:,2));%绘制附件数据图形

%进行灵敏性分析时取消注释
% for de=1.5:0.1:2.5
% for de1=60:70
F=zeros(size(fujian,1),1);
delta=0.5/0.02; %每隔0.5s绘制的间隔
k1=1;k2=1;k3=1;

%对a和k1/k进行大致范围的步长遍历求解
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
            FF(k1,k2)=(fujian(k1,2)+273-xx(k1,k2)).^2;%最小二乘
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

%选择最小二乘最优解
for i=1:size(xxx,2)
    [fm(i),ind(i)]=min(FFF{1,i});
    [fmin,indx]=min(fm);
end

%进行灵敏性分析时取消注释
%         error(d1,d2)=fmin;
%         d1=d1+1;
%     end
%     d1=1;
%     d2=d2+1;
% end

%绘制拟合结果
figure;
plot(fujian(:,1),xxx{indx}(:,ind(indx))-273,'Linewidth',1.2);
hold on;
plot(fujian(:,1),fujian(:,2),'Linewidth',1.2);
grid on;
ylabel('电子板中心区域温度/℃');
xlabel('时间/s');
title('炉温曲线最小二乘拟合曲线图');
legend('拟合结果','附件数据');

%输出参数拟合结果
aaa=0.000015:0.000001:0.000025;
kkk=3000:100:8000;

disp('根据最小二乘法拟合得到的参数a为：');
disp(aaa(indx));
disp('根据最小二乘法拟合得到的参数k1/k为：')
disp(kkk(ind(indx)));

%灵敏性分析图形
% mesh(de,de1,error,'Linewidth',1.2);
% title('温度过渡区域灵敏度分析结果');
% xlabel('冷却区之前的温度过渡区域距离/m','Rotation',0);
% ylabel('冷却区之后的温度过渡区域距离/m');
% zlabel('误差平方和');