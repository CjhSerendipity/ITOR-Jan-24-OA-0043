%%%一种成分血的生产库存分配情况：(以血小板保质期为5天)
clc;
clear all;

tmax=35;%%策略周期
st=2;%%%交叉配型释放期
%%随机产生δ~N（40,15）
u1=50;
s1=10;
Q=normrnd(u1,s1,[1,34]);%%%生产分布
aik=1; %%%每单位产量
q1=Q.*aik;
d1=normrnd(40,15,[1,34]);
%%每天的库存水平
VB1=zeros(1,tmax);
VH1=zeros(1,tmax);
%%库存水平状态
VB11=zeros(1,5);
VH11=zeros(1,5);
cir_1=zeros(tmax,5);
WB1=zeros(1,tmax);
WH1=zeros(1,tmax);
SB1=zeros(1,tmax);
SH1=zeros(1,tmax);
bqh11=zeros(1,5);
%%%初始库存，为1*5的矩阵
%%%%第一天的状态
R1=2.34*s1+u1;
ss1=2.34*s1;
CR=30;
CH=1;
CB=3;
CW=10;
a=0.4;
qh1=sqrt(CR*u1/(2*CH));
VH1(1)=sum(VH11(1:5));
bq11=zeros(1:5);

%%%第一天的情况
D1(1)=0.*(VH1(1)>R1)+qh1.*(ss1<=VH1(1)-R1&VH1(1)-R1<R1)+(qh1+u1)*(VH1(1)<ss1);%%%判断医院想订购量是多少。
VB11(1,5)=VB11(1,5)+q1(1);%%%血液中心第一天的库存水平
VB1(1)=sum(VB11(1:5));
bqh1=D1.*(D1(1)<=VB1(1))+VB1(1).*(VB1(1)<D1(1));%%%确定血液中心分配给医院的分配量是多少。
SB1(1)=(D1(1)-VB1(1)).*(D1(1)>VB1(1));%%%血液中心的短缺量
for i=1:5 %%%%求出分配给医院的量还差多少
sample=sum(VB11(1:i));
temp=min(bqh1,sample);
if temp>=bqh1
    t=bqh1-sum(VB11(1:i-1));
break
end 
end
%%%分配给医院的库龄为多少 量为多少
for p=1:i-1
bq11(1:p)=bq11(1:p)+VB11(1:p);
VB11(1,p)=0;%%%将血液中心已经分配的置零
end
bq11(1:i)=bq11(1:i)+t;
VB11(i)=VB11(i)-t;
VB1(1)=sum(VB11(1:5));%%%血液中心剩余库存水平
WB1(1)=WB1(1)+VB11(1);%%血液中心剩余保质期为0的过期报废。

%%%医院的决策
VH11(1:5)=VH11(1:5)+bq11(1:5);%%%%医院期初库存水平
VH1(1)=sum(VB11(1:5));%%%医院期初库存水平
SH1(1)=SH1(1)+(d1(1)-VH1(1)).*(d1(1)>VH1(1));%%%医院短缺量
%%%计算剩余量
for m=1:5
sample=sum(VH11(1:m));
temp=min(d1(1),sample);
if temp>=d1(1)
    x=d1(1)-sum(VH11(1:m-1));
break
end
end
%%%计算交叉配型备血库存输注完后 的剩余量
for p=1:i-1
cir_1(1,p)=cir_1(1,p)+(1-a)*VH1(1,p);
VH11(1,p)=0;
end
cir_1(1,i)=cir_1(1,i)+(1-a)*x;
VH11(1,i)=VH11(1,i)-x;
WH1(1)=WH1(1)+VH11(1)+cir_1(1);

%%%更新库存状态：
for j=1:4
VB11(1,j)=VB11(1,j+1);
VH11(1,j)=VH11(1,j+1);
cir_1(1,j)=cir_1(1,j+1);
end
VB11(1,5)=0;
VH11(1,5)=0;
cir_1(1,5)=0;
VH1(2)=sum(VH11(1:5));%%%医院期末库存水平
VB1(2)=sum(VB11(1:5));%%%血液中心期末库存水平

%%%从第2天到tmax那天的库存状态
for day=2:tmax
 if day>st
   for j=1:5
      VH11(1,j)=VH11(1,j)+cir_1(day-st,j);
   end
 end
VH1(day)=sum(VH11(1:5));
D1(day)=0.*(VH1(day)>R1)+qh1.*(ss1<=VH1(day)-R1&VH1(day)-R1<R1)+(qh1+u1)*(VH1(day)<ss1);%%%判断医院想订购量是多少。
VB11(1,5)=VB11(1,5)+q1(day);%%%血液中心第t天的库存水平
VB1(day)=sum(VB11(1:5));
bqh1=D1.*(D1(day)<=VB1(day))+VB1(day).*(VB1(day)<D1(day));%%%确定血液中心分配给医院的分配量是多少。
SB1(day)=(D1(day)-VB1(day)).*(D1(day)>VB1(day));%%%血液中心的短缺量

for i=1:5
sample=sum(VB11(1:i));
temp=min(bqh1,sample);
  if temp>=bqh1
     t=bqh1-sum(VB11(1:i-1));
     break
  end
end
for p=1:i-1
bqh11(1:p)=bqh11(1:p)+VB11(1:p);
VB11(1:p)=0;
end
bqh11(1,i)=bqh11(1,i)+t;
VB11(1,i)=VB11(1,i)-t;

WB1(day)=WB1(day)+VB11(1);%%血液中心剩余保质期为0的过期报废。
VB1(day)=sum(VB11(2:5));%%%血液中心剩余库存水平

%%%医院的决策
VH11(1:5)=VH11(1:5)+bqh11(1:5);%%%%医院期初库存水平
VH1(day)=sum(VB11(1:5));%%%医院期初库存水平
SH1(day)=(d1(day)-VH1(day)).*(d1(day)>VH1(day));%%%医院短缺量
for m=1:5
sample=sum(VH11(1:m));
temp=min(d1(day),sample);
if temp>=d1(day)
    x=d1(day)-sum(VH11(1:m-1));
    break
end
end 
%%%计算交叉配型备血库存输注完后 的剩余量
for p=1:i-1
cir_1(day,1:p)=cir_1(day,1:p)+(1-a)*VH1(1,p);
VH11(1,1:p)=0;
end
cir_1(day,i)=cir_1(day,i)+(1-a)*x;
VH11(1,i)=VH11(1,i)-x;
WH1(day)=VH11(1)+cir_1(1);
%%%更新库存状态：
for j=1:4
VB11(1,j)=VB11(1,j+1);
VH11(1,j)=VH11(1,j+1);
  for i=1:tmax
  cir_1(i,j)=cir_1(i,j+1);
  end
end
VB11(1,5)=0;
VH11(1,5)=0;
cir_1(1,5)=0;

end

%%%目标函数
obj_1=(sum(SB1(1:tmax))+sum(SH1(1:tmax)))\sum(q1(1:tmax:1));%%%短缺率最小
obj_2=CB*sum(VB1(1:tmax))+CH*sum(VH1(1:tmax))+CW*(sum(WB1(1:tmax))+sum(WH1(1:tmax)));%%%总成本最低



