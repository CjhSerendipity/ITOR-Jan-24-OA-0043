clear all
clc
tic
global a_kp T X_1t X_2t L_p V_l0Bp V_l0Hp u_tp sigema_tp alpha O_tp gama substitue T_maxtp C_kp C_EBp C_EHp C_Sp C_Wp C_Tp C_ypp
%����
drawing_flag = 1;
data_a_kp = xlsread('data/a_kp.xlsx');
parametes = xlsread('data/parametes.xlsx');
substitue = xlsread('data/substitute.xlsx');
D_tttt = xlsread('data/D_t.xlsx');
for a = 1:7
    for b = 1:7
        if isnan(substitue(a,b))
            substitue(a,b) = inf;
        end
    end
end
Vp = xlsread('data/V.xlsx');
T = 30;
V_l0Bp = zeros(35,7); %ֻ����30�죬���ڳ����ֻ��5��Ŀ��䣬������ѪҺ���Ŀ��伴��35��
V_l0Bp(1:5,1:7) = Vp; %�ڳ�ѪҺ���Ŀ��
V_l0Hp = zeros(35,7); %�ڳ�ҽԺ���
a_kp = data_a_kp(1:11,1:7); %���Բ���
L_p = parametes(1,:); %����
X_1t(1:T) = 200; %tʱȫѪ�ɼ���
X_2t(1:T) = 50; %tʱ������
V_maxBp = inf;%ѪҺ����p�ͳɷ�Ѫ��߿������ 
V_maxHp = inf;%ҽԺp�ͳɷ�Ѫ��߿������ 
T_max = parametes(7,:); 
for t = 1:T
    T_maxtp(t,:) = T_max; %tʱ�ɵ�������
end
u_tp = parametes(8,:); %����
sigema_tp = parametes(9,:); %�����׼��
alpha = parametes(11,:); %��ע����
gama = parametes(10,:); %���������ͷ���
C_kp = data_a_kp(:,8); %���Գɱ�
C_EBp = parametes(2,:); %ѪҺ���Ŀ��ɱ�
C_EHp = parametes(3,:); %ҽԺ���ɱ�
C_Sp = parametes(5,:); %��ȱ�ɱ�
C_Wp = parametes(4,:); %���ϳɱ�
C_Tp = parametes(6,:); %�����ɱ�
C_ypp = substitue; %����ɱ�
O_tp = D_tttt;
nVar=T*11;
VarSize=[1 nVar];
GreyWolves_num=100;
MaxIt=400;  %��������
Archive_size=100;   % �洢��С

alphaa=0.1;  %�������Ͳ���
nGrid=10;   %ÿ��ά�ȵ�������
beta=4;     %�쵼ѡ��ѹ������
gamma=2;    %���⣨��ɾ�����洢���Աѡ��ѹ�� 

GreyWolves=CreateEmptyParticle(GreyWolves_num);


for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    angents_pos.x_kt = initial_population(1);
    angets_obj = objective(angents_pos,1);
    x = angets_obj(1);
    y = angets_obj(2);
    z = angets_obj(3);
    transposition_pos = (angents_pos.x_kt)';
    new_pos(1,:) =  transposition_pos(:);
    GreyWolves(i).Cost = [x;y;z];
    GreyWolves(i).Position = new_pos(1,:);
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination(GreyWolves);
Archive=GetNonDominatedParticles(GreyWolves);
Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alphaa);

for i=1:numel(Archive)
    [Archive(i).GridIndex,Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

for it=1:MaxIt
    a=2-it*((2)/MaxIt);
    popc = GreyWolves;
    for g=1:GreyWolves_num
        
        clear rep2
        clear rep3
        
        %ѡ��alpha,beta,delta����
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        %����ӵ����Ⱥ���еĽ�������������ڶ����ӵ����Ⱥ��Ҳ����ѡ�������쵼�ߡ� 
        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Beta=SelectLeader(rep2,beta);
        end
        
        %����ڶ����ӵ����Ⱥ����һ���⣬��˳�������ͬ�ģ����Ӧ�ӵ������ӵ����Ⱥ����ѡ��delta�쵼�ߡ� 
        if size(Archive,1)>2
            counter=0;
            for newi=1:size(rep2,1)
                if sum(Beta.Position~=rep2(newi).Position)~=0
                    counter=counter+1;
                    rep3(counter,1)=rep2(newi);
                end
            end
            Alpha=SelectLeader(rep3,beta);
        end
        c=2.*rand(1, nVar);
        D=abs(c.*Delta.Position-GreyWolves(g).Position);
        A=2.*a.*rand(1, nVar)-a;
        X1=Delta.Position-A.*abs(D);
        c=2.*rand(1, nVar);
        D=abs(c.*Beta.Position-GreyWolves(g).Position);
        A=2.*a.*rand()-a;
        X2=Beta.Position-A.*abs(D);
        c=2.*rand(1, nVar);
        D=abs(c.*Alpha.Position-GreyWolves(g).Position);
        A=2.*a.*rand()-a;
        X3=Alpha.Position-A.*abs(D);
        GreyWolves(g).Position=(X1+X2+X3)./3;
        aaa = (reshape(GreyWolves(g).Position,11,30))';
        for t = 1:1:T
            x_k = abs(aaa(t,:));
            s_x = sum(x_k(1:4));
            s_x2 = sum(x_k(5:6));
            for h = 1:4
                x_k(h) = fix((x_k(h)/s_x)*X_1t(t));
            end
            for h = 1:2
                if fix(s_x2) >= X_2t(t)
                    x_k(4+h) = fix((x_k(4+h)/s_x2)*X_2t(t));
                else
                    x_k(4+h) = fix(x_k(4+h));
                end  
            end
            s_x3 = sum(x_k(7:8));
            if fix(x_k(7)+x_k(8)) > x_k(2)
                for h = 1:2
                    x_k(6+h) = fix((x_k(6+h)/s_x3)*x_k(2));
                end
            else
                x_k(7) = fix(x_k(7));
                x_k(8) = fix(x_k(8));
            end
            if fix(x_k(9)) > x_k(1)+x_k(8)+x_k(3)
                x_k(9) = x_k(1)+x_k(8)+x_k(3);
            else
                x_k(9) = fix(x_k(9));
            end
            if fix(x_k(10)) > x_k(1)+x_k(2)+x_k(3) + 3*x_k(5)
                x_k(10) = x_k(1)+x_k(2)+x_k(3) + 3*x_k(5);
            end
                x_k(10) = fix(x_k(10));
            if fix(x_k(11)) > x_k(3)
                x_k(11) = x_k(3);
            else
                x_k(11) = fix(x_k(11));
            end
            aaa(t,:) = x_k;
        end
        abc.x_kt = aaa;
        aaaobj = objective(abc,1);
        x = aaaobj(1);
        y = aaaobj(2);
        z = aaaobj(3);
        GreyWolves(g).Cost = [x;y;z];
        transposition_pos = (aaa)';
        GreyWolves(g).Position(1,:) =  transposition_pos(:);
    end
%     %�ϲ���Ⱥ�������¼����֧��ȼ���ӵ����
%     newpop = [GreyWolves; popc];
%     [GreyWolves,F] = nondominatedsort(newpop);
%     GreyWolves = calcrowdingdistance(GreyWolves,F);
%     %����
%     GreyWolves = Sortpop(GreyWolves);
%     %��̭
%     GreyWolves = GreyWolves(1: GreyWolves_num);
%     [GreyWolves,F] = nondominatedsort(GreyWolves);
%     GreyWolves = calcrowdingdistance(GreyWolves,F);
%     GreyWolves = Sortpop(GreyWolves);
%     F1 = GreyWolves(F{1});
%     for f=1:1:size(F1,1)
%         F1obj(f,:) = (F1(f).Cost)';
%     end
%     F1obj = unique(F1obj,'rows');
%     aa = size(F,2);
%     disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1)) '��֧��ȼ� ' num2str(aa)]);
%     % ��ͼ
%     plotpp_Chinese(F1obj);
%     pause(0.01);

    
    GreyWolves=DetermineDomination(GreyWolves);
    non_dominated_wolves=GetNonDominatedParticles(GreyWolves);
    
    Archive=[Archive
        non_dominated_wolves];
    
    Archive=DetermineDomination(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex,Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alphaa);
        
    end
    
%     [Grey,F] = nondominatedsort(GreyWolves);
%     AAA = size(F,2);
    for f=1:1:size(Archive,1)
        F1obj(f,:) = (Archive(f).Cost)';
    end
    F1obj = unique(F1obj,'rows');
    %disp(['Iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(size(F1obj,1)) '��֧��ȼ� ' num2str(AAA)]);
    disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(size(F1obj,1))]);
    %aa = size(Archive,1);
    % ��ͼ
    plotpp_Chinese(F1obj);
    pause(0.01);
    it = it+1;

end
toc
% for i =1:size(Archive,1)
%     pop(i).value = Archive(i).Position;
%     pop(i).obj = Archive(i).Cost;
% end
% for i =1:size(GreyWolves,1)
%     popc(i).value = GreyWolves(i).Position;
%     popc(i).obj = GreyWolves(i).Cost;
% end
% [popc,f] = nondominatedsort(popc);
%[pop,F] = nondominatedsort(pop);
% if size(F,1)
%     hv = HV(F1obj );
% end
disp(['����ʱ��: ',num2str(toc)]);

time = toc;


