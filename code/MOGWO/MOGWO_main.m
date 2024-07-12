clear all
clc
tic
global a_kp T X_1t X_2t L_p V_l0Bp V_l0Hp u_tp sigema_tp alpha O_tp gama substitue T_maxtp C_kp C_EBp C_EHp C_Sp C_Wp C_Tp C_ypp
%参数
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
V_l0Bp = zeros(35,7); %只计算30天，且期初最大只有5天的库龄，故所有血液最大的库龄即到35天
V_l0Bp(1:5,1:7) = Vp; %期初血液中心库存
V_l0Hp = zeros(35,7); %期初医院库存
a_kp = data_a_kp(1:11,1:7); %策略产量
L_p = parametes(1,:); %库龄
X_1t(1:T) = 200; %t时全血采集量
X_2t(1:T) = 50; %t时单采量
V_maxBp = inf;%血液中心p型成分血最高库存限制 
V_maxHp = inf;%医院p型成分血最高库存限制 
T_max = parametes(7,:); 
for t = 1:T
    T_maxtp(t,:) = T_max; %t时可调剂上限
end
u_tp = parametes(8,:); %需求
sigema_tp = parametes(9,:); %需求标准差
alpha = parametes(11,:); %输注比例
gama = parametes(10,:); %交叉配型释放期
C_kp = data_a_kp(:,8); %策略成本
C_EBp = parametes(2,:); %血液中心库存成本
C_EHp = parametes(3,:); %医院库存成本
C_Sp = parametes(5,:); %短缺成本
C_Wp = parametes(4,:); %报废成本
C_Tp = parametes(6,:); %调剂成本
C_ypp = substitue; %替代成本
O_tp = D_tttt;
nVar=T*11;
VarSize=[1 nVar];
GreyWolves_num=100;
MaxIt=400;  %迭代次数
Archive_size=100;   % 存储大小

alphaa=0.1;  %网格膨胀参数
nGrid=10;   %每个维度的网格数
beta=4;     %领导选择压力参数
gamma=2;    %额外（待删除）存储库成员选择压力 

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
        
        %选择alpha,beta,delta灰狼
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        %如果最不拥挤的群体中的解少于三个，则第二个最不拥挤的群体也可以选择其他领导者。 
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
        
        %如果第二个最不拥挤的群体有一个解，则此场景是相同的，因此应从第三个最不拥挤的群体中选择delta领导者。 
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
%     %合并种群，并重新计算非支配等级和拥挤度
%     newpop = [GreyWolves; popc];
%     [GreyWolves,F] = nondominatedsort(newpop);
%     GreyWolves = calcrowdingdistance(GreyWolves,F);
%     %排序
%     GreyWolves = Sortpop(GreyWolves);
%     %淘汰
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
%     disp(['Iteration ' num2str(it) ': Number of F1 Members = ' num2str(numel(F1)) '非支配等级 ' num2str(aa)]);
%     % 绘图
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
    %disp(['Iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(size(F1obj,1)) '非支配等级 ' num2str(AAA)]);
    disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(size(F1obj,1))]);
    %aa = size(Archive,1);
    % 绘图
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
disp(['运行时间: ',num2str(toc)]);

time = toc;


