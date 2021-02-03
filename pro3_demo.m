clc
clear
%绘制问题三结果图像

%参数设置
aa=0.000016;
kk=7000;
de=2/100;
de1=66/100;
v=8.59/6/100;
dt=0.02;

%温度设置
T0=25;
T1=181.5906;
T2=193.1463;
T3=226.6849;
T4=264.3210;

%速度设置
v=85.9851/60/100;

%温度模型
Tair=[];
S_before=0.25;
S_after=11*0.305+10*0.05+0.25;
t=-ceil(S_before/v*50)/50:0.02:ceil(S_after/v*50)/50;
for i=1:size(t,2)
    if(v*t(i)<-de)
        Tair(i)=T0+273;
    elseif(v*t(i)>=-de&&v*t(i)<=de)
        Tair(i)=(v*t(i)+de).*(T1-T0)./(2.*de)+T0+273;
    elseif(v*t(i)>de&&v*t(i)<=1.725-de)
        Tair(i)=T1+273;
    elseif(v*t(i)>1.725-de&&v*t(i)<=1.775+de)
        Tair(i)=(v*t(i)+de-1.725).*(T2-T1)./(2.*de+0.05)+T1+273;
    elseif(v*t(i)>1.775+de&&v*t(i)<=2.08-de)
        Tair(i)=T2+273;
    elseif(v*t(i)>2.08-de&&v*t(i)<=2.13+de)
        Tair(i)=(v*t(i)+de-2.08).*(T3-T2)./(2.*de+0.05)+T2+273;
    elseif(v*t(i)>2.13+de&&v*t(i)<=2.435-de)
        Tair(i)=T3+273;
    elseif(v*t(i)>2.435-de&&v*t(i)<=2.485+de)
        Tair(i)=(v*t(i)+de-2.435).*(T4-T3)./(2.*de+0.05)+T3+273;
    elseif(v*t(i)>2.485+de&&v*t(i)<=3.145-de)
        Tair(i)=T4+273;
    elseif(v*t(i)>3.145-de&&v*t(i)<=3.195+de1)
        Tair(i)=(v*t(i)+de-3.145).*(T0-T4)./(de+de1+0.05)+T4+273;
    else
        Tair(i)=T0+273;
    end
end
Tair=Tair';

%空间步长
dx=0.005/1000;
d=0.15/1000;
nn=round(d/dx);
T=zeros(size(Tair,1),round(d/dx)+1);
T(1,:)=25+273;

%求解差分方程
for n=1:size(Tair,1)-1
    for i=2:nn
        T(n+1,i)=(T(n,i+1)-2.*T(n,i)+T(n,i-1)).*dt.*(aa.^2)./(dx.^2)+T(n,i);
    end
    T(n+1,1)=(dx.*kk.*Tair(n+1,1)+T(n+1,2))./(kk.*dx+1);
    T(n+1,nn+1)=(dx.*kk.*Tair(n+1,1)+T(n+1,nn))./(kk.*dx+1);
end

%从30℃开始绘制
t=t-t(1,1);
for i=1:size(T,1)
    if(T(i,16)>30+273)
        n=i;
        break;
    end
end

%每隔0.5S绘制
delta=0.5/0.02;
k1=1;
for i=n:delta:size(T,1)
    xx(k1,1)=T(i,16);
    tt(k1,1)=t(i);
    k1=k1+1;
end
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
plot([0,tt(endind)],[217,217],'k--','Linewidth',1.2);
plot([tt(endind),tt(endind)],[0,endmax-273],'k--','Linewidth',1.2);
text(3,endmax-273+8,'最高峰值','FontSize',10);
text(3,217-8,'217℃','FontSize',10);
xx=xx-273;
for i=1:size(t217:endind,2)
    hold on;
    plot([tt(t217+i) tt(t217+i)],[217 xx(t217+i)],'k:');
end
