clc
clear
%问题一求解

%各个参数设置
aa=0.000016;
kk=7000;
de=2/100;
de1=66/100;
v=7.8/6/100;
dt=0.02;

%温度模型
Tair=[];
S_before=0.25;
S_after=11*0.305+10*0.05+0.25;
t=-ceil(S_before/v*50)/50:0.02:ceil(S_after/v*50)/50;
for i=1:size(t,2)
    if(v*t(i)<-de)
        Tair(i)=25+273;
    elseif(v*t(i)>=-de&&v*t(i)<=de)
        Tair(i)=(v*t(i)+de).*(173-25)./(2.*de)+25+273;
    elseif(v*t(i)>de&&v*t(i)<=1.725-de)
        Tair(i)=173+273;
    elseif(v*t(i)>1.725-de&&v*t(i)<=1.775+de)
        Tair(i)=(v*t(i)+de-1.725).*(198-173)./(2.*de+0.05)+173+273;
    elseif(v*t(i)>1.775+de&&v*t(i)<=2.08-de)
        Tair(i)=198+273;
    elseif(v*t(i)>2.08-de&&v*t(i)<=2.13+de)
        Tair(i)=(v*t(i)+de-2.08).*(230-198)./(2.*de+0.05)+198+273;
    elseif(v*t(i)>2.13+de&&v*t(i)<=2.435-de)
        Tair(i)=230+273;
    elseif(v*t(i)>2.435-de&&v*t(i)<=2.485+de)
        Tair(i)=(v*t(i)+de-2.435).*(257-230)./(2.*de+0.05)+230+273;
    elseif(v*t(i)>2.485+de&&v*t(i)<=3.145-de)
        Tair(i)=257+273;
    elseif(v*t(i)>3.145-de&&v*t(i)<=3.195+de1)
        Tair(i)=(v*t(i)+de-3.145).*(25-257)./(de+de1+0.05)+257+273;
    else
        Tair(i)=25+273;
    end
end
Tair=Tair';

%步长设置
dx=0.005/1000;
d=0.15/1000;
nn=round(d/dx);
T=zeros(size(Tair,1),round(d/dx)+1);
T(1,:)=25+273;

%有限差分法求解热传导模型
for n=1:size(Tair,1)-1
    for i=2:nn
        T(n+1,i)=(T(n,i+1)-2.*T(n,i)+T(n,i-1)).*dt.*(aa.^2)./(dx.^2)+T(n,i);
    end
    T(n+1,1)=(dx.*kk.*Tair(n+1,1)+T(n+1,2))./(kk.*dx+1);
    T(n+1,nn+1)=(dx.*kk.*Tair(n+1,1)+T(n+1,nn))./(kk.*dx+1);
end

%从升温至30℃以上的点开始绘制
t=t-t(1,1);
for i=1:size(T,1)
    if(T(i,16)>30+273)
        n=i;
        break;
    end
end
delta=0.5/0.02;
k1=1;

%每隔0.5S进行绘制
for i=n:delta:size(T,1)
    xx(k1,1)=T(i,16);
    tt(k1,1)=t(i);
    k1=k1+1;
end

%绘制结果
plot(tt,xx-273,'-','Linewidth',1.5);
axis([0,350,280-273,550-273]);
grid on;
title('炉温曲线');
xlabel('时间/s');
ylabel('电子板中心区域温度/℃');
[endmax,endind]=max(xx);
hold on;
plot([0,tt(endind)],[endmax-273,endmax-273],'k--','Linewidth',1.2);
for i=1:size(xx,1)
    if(xx(i)>217+273)
        t217=i;
        break;
    end
end
plot([0,tt(t217)],[217,217],'k--','Linewidth',1.2);
text(3,endmax-273+8,'最高峰值','FontSize',10);
text(3,217-8,'217℃','FontSize',10);
xx=xx-273;


%输出指标
disp('小温区3中点焊接区域中心的温度：')
disp(T(4279,16)-273);
disp('小温区6中点焊接区域中心的温度：')
disp(T(8375,16)-273);
disp('小温区7中点焊接区域中心的温度：')
disp(T(9740,16)-273);
disp('小温区8结束处焊接区域中心的温度：')
disp(T(11692,16)-273);