clc
clear all;
tic
global a_kp T X_1t X_2t L_p V_l0Bp V_l0Hp u_tp sigema_tp alpha O_tp gama substitue T_maxtp C_kp C_EBp C_EHp C_Sp C_Wp C_Tp C_ypp
%参数
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
% D_tHp = zeros(30,7);
% for ii = 1:7
%     D_tHp(:,ii) = abs(fix(normrnd(u_tp(ii),sigema_tp(ii),[30,1])));%t天初需求量
% end
O_tp = D_tttt;
%算法相关参数
SearchAgents_no = 100;
Max_iter = 400;
%生成初始种群，计算目标值
angents_pos = initial_population(SearchAgents_no);
angents_obj = Objective(angents_pos,SearchAgents_no);
% initialize position vector and score for the leader
Leader_pos = zeros(1,T*size(a_kp,1));
Leader_score = inf; 
Convergence_curve = zeros(1,Max_iter);
tt  = 0;% Loop counter
%非支配排序准备工作
empty.value = [];
empty.obj = [];
empty.rank = [];
empty.domination = [];
empty.dominated = 0;
empty.crowdingdistance = [];
pop = repmat(empty, SearchAgents_no, 1);
for p = 1:1:SearchAgents_no
    x = angents_obj(p).obj1;
    y = angents_obj(p).obj2;
    z = angents_obj(p).obj3;
    transposition_pos = (angents_pos(p).x_kt)';
    new_pos(p,:) =  transposition_pos(:);
    pop(p).obj = [x;y;z];
    pop(p).value = new_pos(p,:);
end
%非支配排序
[pop,F] = nondominatedsort(pop);
%拥挤度计算
pop = calcrowdingdistance(pop,F);
%主程序
while tt < Max_iter
    pop = Sortpop(pop);
    popc = repmat(empty, SearchAgents_no ,1);
%     for i = 2:size(SearchAgents_no,1)
%         angents_pos(i-1,:) = pop(i).value;
%         angents_obj(i-1,:) = pop(i).obj;
%     end
    a = 2-t*((2)/Max_iter); 
    a2 = -1+t*((-1)/Max_iter);
    F1 = pop(F{1});
    % Update the Position of search agents 
    for i = 1:size(new_pos,1)
        f = randi(numel(F1));
        leader_pop = F1(f);
        Leader_pos = leader_pop.value;
        Leader_score = leader_pop.obj;
        r1 = rand(); 
        r2 = rand(); 
        A = 2*a*r1-a;  
        C = 2*r2;      
        b = 1;               
        l = (a2-1)*rand+1;   
        p = rand();    
        for j = 1:size(new_pos,2)
            if p < 0.5   
                if abs(A) >= 1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = new_pos(rand_leader_index, :);
                    D_X_rand = abs(C*X_rand(j) - new_pos(i,j)); 
                    new_pos(i,j) = X_rand(j)-A*D_X_rand;    
                elseif abs(A) < 1
                    D_Leader = abs(C*Leader_pos(j) - new_pos(i,j)); 
                    new_pos(i,j) = Leader_pos(j)-A*D_Leader;     
                end
            elseif p>=0.5
                distance2Leader = abs(Leader_pos(j) - new_pos(i,j));
                new_pos(i,j) = distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
            end
        end
        angents_pos(i).x_kt = (reshape(new_pos(i,:),11,30))';
        for t = 1:1:T
            x_k = abs(angents_pos(i).x_kt(t,:));
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
            angents_pos(i).x_kt(t,:) = x_k;
        end   
    end
    angents_obj = Objective(angents_pos,SearchAgents_no);
    for p = 1:1:SearchAgents_no
        x = angents_obj(p).obj1;
        y = angents_obj(p).obj2;
        z = angents_obj(p).obj3;
        transposition_pos = (angents_pos(p).x_kt)';
        new_pos(p,:) =  transposition_pos(:);
        popc(p).obj = [x;y;z];
        popc(p).value = new_pos(p,:);
    end
    %合并种群，并重新计算非支配等级和拥挤度
    newpop = [pop; popc];
    [pop,F] = nondominatedsort(newpop);
    pop = calcrowdingdistance(pop,F);
    %排序
    pop = Sortpop(pop);
    %淘汰
    pop = pop(1: SearchAgents_no);
    [pop,F] = nondominatedsort(pop);
    pop = calcrowdingdistance(pop,F);
    pop = Sortpop(pop);
    F1 = pop(F{1});
    for f=1:1:size(F1,1)
        F1obj(f,:) = (F1(f).obj)';
    end
    F1obj = unique(F1obj,'rows');
    aa = size(F,2);
    disp(['Iteration ' num2str(tt) ': Number of F1 Members = ' num2str(numel(F1)) '非支配等级 ' num2str(aa)]);
    % 绘图
    plotpp_Chinese(F1);
    pause(0.01);
    tt = tt+1;
end
hv = HV(F1obj);                    
toc
disp(['运行时间: ',num2str(toc)]);
time = toc;
