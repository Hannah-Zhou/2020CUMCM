clc
clear
%问题二求解

%参数设置
aa=0.000016;
kk=7000;
de=2/100;
de1=66/100;

k=1;
%粗遍历和精遍历
for v_=65:0.2:100  %粗遍历
    %for v_=72.4:0.01:74.4  %精遍历
    v=v_/6/1000;
    dt=0.02;%时间步长
    
    %温度模型
    Tair=[];
    S_before=0.25;
    S_after=11*0.305+10*0.05+0.25;
    t=-ceil(S_before/v*50)/50:0.02:ceil(S_after/v*50)/50;
    for i=1:size(t,2)
        if(v*t(i)<-de)
            Tair(i)=25+273;
        elseif(v*t(i)>=-de&&v*t(i)<=de)
            Tair(i)=(v*t(i)+de).*(182-25)./(2.*de)+25+273;
        elseif(v*t(i)>de&&v*t(i)<=1.725-de)
            Tair(i)=182+273;
        elseif(v*t(i)>1.725-de&&v*t(i)<=1.775+de)
            Tair(i)=(v*t(i)+de-1.725).*(203-182)./(2.*de+0.05)+182+273;
        elseif(v*t(i)>1.775+de&&v*t(i)<=2.08-de)
            Tair(i)=203+273;
        elseif(v*t(i)>2.08-de&&v*t(i)<=2.13+de)
            Tair(i)=(v*t(i)+de-2.08).*(237-203)./(2.*de+0.05)+203+273;
        elseif(v*t(i)>2.13+de&&v*t(i)<=2.435-de)
            Tair(i)=237+273;
        elseif(v*t(i)>2.435-de&&v*t(i)<=2.485+de)
            Tair(i)=(v*t(i)+de-2.435).*(254-237)./(2.*de+0.05)+237+273;
        elseif(v*t(i)>2.485+de&&v*t(i)<=3.145-de)
            Tair(i)=254+273;
        elseif(v*t(i)>3.145-de&&v*t(i)<=3.195+de1)
            Tair(i)=(v*t(i)+de-3.145).*(25-254)./(de+de1+0.05)+254+273;
        else
            Tair(i)=25+273;
        end
    end
    Tair=Tair';
    
    %空间步长
    dx=0.005/1000;
    d=0.15/1000;
    nn=round(d/dx);
    T=zeros(size(Tair,1),round(d/dx)+1);
    T(1,:)=25+273;
    
    
    %有限差分法求解差分方程
    for n=1:size(Tair,1)-1
        for i=2:nn
            T(n+1,i)=(T(n,i+1)-2.*T(n,i)+T(n,i-1)).*dt.*(aa.^2)./(dx.^2)+T(n,i);
        end
        T(n+1,1)=(dx.*kk.*Tair(n+1,1)+T(n+1,2))./(kk.*dx+1);
        T(n+1,nn+1)=(dx.*kk.*Tair(n+1,1)+T(n+1,nn))./(kk.*dx+1);
    end
    
    %制程界限计算
    t=t-t(1,1);
    
    %斜率 用差分代替
    for i=1:size(T,1)-1
        Tv(i)=(T(i+1,16)-T(i,16))/dt;
        
    end
    %150℃
    for i=1:size(T,1)
        if T(i,16)>=150+273
            t1=i;
            break;
        end
    end
    %190℃
    for i=1:size(T,1)
        if T(i,16)>=190+273
            t2=i;
            break;
        end
    end
    deltat=(t2-t1)*dt;
    t3=length(find(T(:,16)>217+273))*dt;
    %峰值
    Tmax=max(T(:,16));
    
    %记录满足制程界限的速度
    if(max(abs(Tv))<=3&&deltat>=60&&deltat<=120&&t3>=40&&t3<=90&&(Tmax>=240+273)&&(Tmax<=250+273))
        v_sel(k)=v*1000*6;
        k=k+1;
    end
end
disp('本次遍历得到的最大速度为');
disp(v_sel(length(v_sel)));