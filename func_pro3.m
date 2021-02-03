%问题三优化目标函数
function f=func_pro3(xx)
%f代表面积

%参数设定
T1=xx(1);
T2=xx(2);
T3=xx(3);
T4=xx(4);
v=xx(5)/60/100;
aa=0.000016;
kk=7000;
de=2/100;
de1=66/100;
dt=0.02;

%温度模型
Tair=[];
S_before=0.25;
S_after=11*0.305+10*0.05+0.25;
t=-ceil(S_before/v*50)/50:0.02:ceil(S_after/v*50)/50;
T0=25;
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

%有限差分法求解热传导方程
for n=1:size(Tair,1)-1
    for i=2:nn
        T(n+1,i)=(T(n,i+1)-2.*T(n,i)+T(n,i-1)).*dt.*(aa.^2)./(dx.^2)+T(n,i);
    end
    T(n+1,1)=(dx.*kk.*Tair(n+1,1)+T(n+1,2))./(kk.*dx+1);
    T(n+1,nn+1)=(dx.*kk.*Tair(n+1,1)+T(n+1,nn))./(kk.*dx+1);
end

%制程界限条件
for i=1:size(T,1)-1
    Tv(i)=(T(i+1,16)-T(i,16))/dt;
end

for i=1:size(T,1)
    if T(i,16)>=150+273
        t1=i;
        break;
    end
end

for i=1:size(T,1)
    if T(i,16)>=190+273
        t2=i;
        break;
    end
end

deltat=(t2-t1)*dt;
t3=length(find(T(:,16)>217+273))*dt;
[Tmax,index]=max(T(:,16));

%求阴影区域面积 用梯形面积代替积分函数 加快求解速度
for i=1:size(T,1)
    if T(i,16)>=217+273
        t4=i;
        break;
    end
end
f=T(t4:index,16)-217-273;
f=0.01.*(2.*sum(f)-f(1)-f(length(f)));

%若不能满足制程界限 就面积赋值为无穷大
if(max(abs(Tv))>3||deltat<60||deltat>120||t3<40||t3>90||(Tmax<240+273)||(Tmax>250+273))
    f=1000000000;
end
end