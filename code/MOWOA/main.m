clc
clear all;
tic
global a_kp T X_1t X_2t L_p V_l0Bp V_l0Hp u_tp sigema_tp alpha O_tp gama substitue T_maxtp C_kp C_EBp C_EHp C_Sp C_Wp C_Tp C_ypp
%����
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
% D_tHp = zeros(30,7);
% for ii = 1:7
%     D_tHp(:,ii) = abs(fix(normrnd(u_tp(ii),sigema_tp(ii),[30,1])));%t���������
% end
O_tp = D_tttt;
%�㷨��ز���
SearchAgents_no = 100;
Max_iter = 400;
%���ɳ�ʼ��Ⱥ������Ŀ��ֵ
angents_pos = initial_population(SearchAgents_no);
angents_obj = Objective(angents_pos,SearchAgents_no);
% initialize position vector and score for the leader
Leader_pos = zeros(1,T*size(a_kp,1));
Leader_score = inf; 
Convergence_curve = zeros(1,Max_iter);
tt  = 0;% Loop counter
%��֧������׼������
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
%��֧������
[pop,F] = nondominatedsort(pop);
%ӵ���ȼ���
pop = calcrowdingdistance(pop,F);
%������
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
    %�ϲ���Ⱥ�������¼����֧��ȼ���ӵ����
    newpop = [pop; popc];
    [pop,F] = nondominatedsort(newpop);
    pop = calcrowdingdistance(pop,F);
    %����
    pop = Sortpop(pop);
    %��̭
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
    disp(['Iteration ' num2str(tt) ': Number of F1 Members = ' num2str(numel(F1)) '��֧��ȼ� ' num2str(aa)]);
    % ��ͼ
    plotpp_Chinese(F1);
    pause(0.01);
    tt = tt+1;
end
hv = HV(F1obj);                    
toc
disp(['����ʱ��: ',num2str(toc)]);
time = toc;
