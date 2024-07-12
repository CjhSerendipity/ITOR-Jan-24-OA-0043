function [angents_obj] = new_Objective(angents_pos,SearchAgents_no)
    global T L_p a_kp V_l0Bp V_l0Hp alpha D_tHp
    P = size(a_kp,2);
    for pp = 1:1:P
        if L_p(pp) > 35
            L_p(pp) = 35;
        end
    end
    for a = 1:1:SearchAgents_no
        obj1 = 0;
        obj2 = 0;
        obj3 = 0;
        x_kt = angents_pos(a).x_kt;
        V_ltBp = zeros(35,P,T); 
        V_ltBp(:,:,1) = V_l0Bp;
        V_ltHp = zeros(35,P,T); 
        V_ltHp(:,:,1) = V_l0Hp;
        q_ltp = zeros(35,P,T+1); %t天初分配量
        W_tBp = zeros(T,P); %血液中心t天末报废量
        W_tHp = zeros(T,P); %医院t天末报废量
        B_ltp = zeros(35,P,T+1); %t天初交叉配型返回量
        S_tp = zeros(T,P); %t天末短缺量
        for t = 1:1:T
            Q_tp(t,:) =  x_kt(t,:)*a_kp; %formula4.1
        end
        for p = 1:1:P
            for t = 1:1:T
                for l = L_p(p):-1:1
                    %formula4.2
                    if (D_tHp(t,p) - sum(V_ltBp(l:L_p(p),p,t))) > 0
                        q_ltp(l,p,t) = V_ltBp(l,p,t);
                    else
                        if  V_ltBp(l,p,t) - (sum(V_ltBp(l:L_p(p),p,t)) - D_tHp(t,p)) > 0
                            q_ltp(l,p,t) = V_ltBp(l,p,t) - (sum(V_ltBp(l:L_p(p),p,t)) - D_tHp(t,p));
                        else
                            q_ltp(l,p,t) = 0;
                        end
                    end
                    if l == 1
                        V_ltBp(l,p,t+1) = Q_tp(t,p);
                        V_ltHp(l,p,t+1) = q_ltp(l,p,t);
                    else
                        V_ltBp(l,p,t+1) = V_ltBp(l-1,p,t) - q_ltp(l-1,p,t);
                        if alpha(p)*D_tHp(t,p) - sum(V_ltHp(l:L_p(p),p,t)) <=0 
                            V_ltHp(l,p,t+1) =  V_ltHp(l-1,p,t) + B_ltp(l,p,t+1) + q_ltp(l,p,t+1);
                        else
                            if V_ltHp(l-1,p,t) - (alpha(p)*D_tHp(t,p) - sum(V_ltHp(l:L_p(p),p,t))) > 0
                                V_ltHp(l,p,t+1) = V_ltHp(l-1,p,t) - (alpha(p)*D_tHp(t,p) - sum(V_ltHp(l:L_p(p),p,t))) + B_ltp(l,p,t+1) + q_ltp(l,p,t+1);
                            else
                                V_ltHp(l,p,t+1) = B_ltp(l,p,t+1) + q_ltp(l,p,t+1);
                            end
                        end
                    end
                end
                if L_p(p) <= 35
                    if V_ltBp(L_p(p),p,t) > q_ltp(L_p(p),p,t)
                        W_tBp(t,p) =  V_ltBp(L_p(p),p,t) - q_ltp(L_p(p),p,t);
                    else
                        W_tBp(t,p) = 0;
                    end
                    if V_ltHp(L_p(p),p,t) > alpha(p)*D_tHp(t,p)
                        W_tHp(t,p) = V_lHBp(L_p(p),p,t) - alpha(p)*D_tHp(t,p);
                    else
                        W_tHp(t,p) = 0;
                    end
                    if D_tHp(t,p) > sum(V_ltHp(:,p,t)) + sum(V_ltBp(:,p,t))
                        S_tp(t,p) = D_tHp(t,p) - (sum(V_ltBp(:,p,t)));
                    else
                        S_tp(t,p) = 0;
                    end
                end
            end
        end
        %调剂或者替代策略选择
        
        obj1 = sum(S_tp,'all')/sum(D_tHp,'all');
        angents_obj(a).obj1 = obj1;
        angents_obj(a).obj2 = obj2;
        angents_obj(a).obj3 = obj3;
    end
end



