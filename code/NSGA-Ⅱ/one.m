%%%һ�ֳɷ�Ѫ�����������������(��ѪС�屣����Ϊ5��)
clc;
clear all;

tmax=35;%%��������
st=2;%%%���������ͷ���
%%���������~N��40,15��
u1=50;
s1=10;
Q=normrnd(u1,s1,[1,34]);%%%�����ֲ�
aik=1; %%%ÿ��λ����
q1=Q.*aik;
d1=normrnd(40,15,[1,34]);
%%ÿ��Ŀ��ˮƽ
VB1=zeros(1,tmax);
VH1=zeros(1,tmax);
%%���ˮƽ״̬
VB11=zeros(1,5);
VH11=zeros(1,5);
cir_1=zeros(tmax,5);
WB1=zeros(1,tmax);
WH1=zeros(1,tmax);
SB1=zeros(1,tmax);
SH1=zeros(1,tmax);
bqh11=zeros(1,5);
%%%��ʼ��棬Ϊ1*5�ľ���
%%%%��һ���״̬
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

%%%��һ������
D1(1)=0.*(VH1(1)>R1)+qh1.*(ss1<=VH1(1)-R1&VH1(1)-R1<R1)+(qh1+u1)*(VH1(1)<ss1);%%%�ж�ҽԺ�붩�����Ƕ��١�
VB11(1,5)=VB11(1,5)+q1(1);%%%ѪҺ���ĵ�һ��Ŀ��ˮƽ
VB1(1)=sum(VB11(1:5));
bqh1=D1.*(D1(1)<=VB1(1))+VB1(1).*(VB1(1)<D1(1));%%%ȷ��ѪҺ���ķ����ҽԺ�ķ������Ƕ��١�
SB1(1)=(D1(1)-VB1(1)).*(D1(1)>VB1(1));%%%ѪҺ���ĵĶ�ȱ��
for i=1:5 %%%%��������ҽԺ�����������
sample=sum(VB11(1:i));
temp=min(bqh1,sample);
if temp>=bqh1
    t=bqh1-sum(VB11(1:i-1));
break
end 
end
%%%�����ҽԺ�Ŀ���Ϊ���� ��Ϊ����
for p=1:i-1
bq11(1:p)=bq11(1:p)+VB11(1:p);
VB11(1,p)=0;%%%��ѪҺ�����Ѿ����������
end
bq11(1:i)=bq11(1:i)+t;
VB11(i)=VB11(i)-t;
VB1(1)=sum(VB11(1:5));%%%ѪҺ����ʣ����ˮƽ
WB1(1)=WB1(1)+VB11(1);%%ѪҺ����ʣ�ౣ����Ϊ0�Ĺ��ڱ��ϡ�

%%%ҽԺ�ľ���
VH11(1:5)=VH11(1:5)+bq11(1:5);%%%%ҽԺ�ڳ����ˮƽ
VH1(1)=sum(VB11(1:5));%%%ҽԺ�ڳ����ˮƽ
SH1(1)=SH1(1)+(d1(1)-VH1(1)).*(d1(1)>VH1(1));%%%ҽԺ��ȱ��
%%%����ʣ����
for m=1:5
sample=sum(VH11(1:m));
temp=min(d1(1),sample);
if temp>=d1(1)
    x=d1(1)-sum(VH11(1:m-1));
break
end
end
%%%���㽻�����ͱ�Ѫ�����ע��� ��ʣ����
for p=1:i-1
cir_1(1,p)=cir_1(1,p)+(1-a)*VH1(1,p);
VH11(1,p)=0;
end
cir_1(1,i)=cir_1(1,i)+(1-a)*x;
VH11(1,i)=VH11(1,i)-x;
WH1(1)=WH1(1)+VH11(1)+cir_1(1);

%%%���¿��״̬��
for j=1:4
VB11(1,j)=VB11(1,j+1);
VH11(1,j)=VH11(1,j+1);
cir_1(1,j)=cir_1(1,j+1);
end
VB11(1,5)=0;
VH11(1,5)=0;
cir_1(1,5)=0;
VH1(2)=sum(VH11(1:5));%%%ҽԺ��ĩ���ˮƽ
VB1(2)=sum(VB11(1:5));%%%ѪҺ������ĩ���ˮƽ

%%%�ӵ�2�쵽tmax����Ŀ��״̬
for day=2:tmax
 if day>st
   for j=1:5
      VH11(1,j)=VH11(1,j)+cir_1(day-st,j);
   end
 end
VH1(day)=sum(VH11(1:5));
D1(day)=0.*(VH1(day)>R1)+qh1.*(ss1<=VH1(day)-R1&VH1(day)-R1<R1)+(qh1+u1)*(VH1(day)<ss1);%%%�ж�ҽԺ�붩�����Ƕ��١�
VB11(1,5)=VB11(1,5)+q1(day);%%%ѪҺ���ĵ�t��Ŀ��ˮƽ
VB1(day)=sum(VB11(1:5));
bqh1=D1.*(D1(day)<=VB1(day))+VB1(day).*(VB1(day)<D1(day));%%%ȷ��ѪҺ���ķ����ҽԺ�ķ������Ƕ��١�
SB1(day)=(D1(day)-VB1(day)).*(D1(day)>VB1(day));%%%ѪҺ���ĵĶ�ȱ��

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

WB1(day)=WB1(day)+VB11(1);%%ѪҺ����ʣ�ౣ����Ϊ0�Ĺ��ڱ��ϡ�
VB1(day)=sum(VB11(2:5));%%%ѪҺ����ʣ����ˮƽ

%%%ҽԺ�ľ���
VH11(1:5)=VH11(1:5)+bqh11(1:5);%%%%ҽԺ�ڳ����ˮƽ
VH1(day)=sum(VB11(1:5));%%%ҽԺ�ڳ����ˮƽ
SH1(day)=(d1(day)-VH1(day)).*(d1(day)>VH1(day));%%%ҽԺ��ȱ��
for m=1:5
sample=sum(VH11(1:m));
temp=min(d1(day),sample);
if temp>=d1(day)
    x=d1(day)-sum(VH11(1:m-1));
    break
end
end 
%%%���㽻�����ͱ�Ѫ�����ע��� ��ʣ����
for p=1:i-1
cir_1(day,1:p)=cir_1(day,1:p)+(1-a)*VH1(1,p);
VH11(1,1:p)=0;
end
cir_1(day,i)=cir_1(day,i)+(1-a)*x;
VH11(1,i)=VH11(1,i)-x;
WH1(day)=VH11(1)+cir_1(1);
%%%���¿��״̬��
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

%%%Ŀ�꺯��
obj_1=(sum(SB1(1:tmax))+sum(SH1(1:tmax)))\sum(q1(1:tmax:1));%%%��ȱ����С
obj_2=CB*sum(VB1(1:tmax))+CH*sum(VH1(1:tmax))+CW*(sum(WB1(1:tmax))+sum(WH1(1:tmax)));%%%�ܳɱ����



