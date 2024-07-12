function [angents_obj] = Objective(angents_pos,SearchAgents_no) 
    global T L_p a_kp V_l0Bp V_l0Hp alpha O_tp gama T_maxtp substitue C_kp C_EBp C_EHp C_Sp C_Wp C_Tp C_ypp
    P = size(a_kp,2);
    for pp = 1:1:P
        if L_p(pp) > 35
            L_p(pp) = 35;
        end
    end
    %for a = 1:1:1
    for a = 1:1:SearchAgents_no
        obj1 = 0;
        obj2 = 0;
        obj3 = 0;
        Tt_maxtp = T_maxtp;
        x_kt = angents_pos(a).x_kt;
        V_ltBp = zeros(35,P,T); 
        V_ltBp(:,:,1) = V_l0Bp; %血液中心期初库存
        V_ltHp = zeros(35,P,T); 
        V_ltHp(:,:,1) = V_l0Hp; %医院期初库存
        q_ltp = zeros(35,P,T+1); %t天初分配量
        W_tBp = zeros(T,P); %血液中心t天末报废量
        W_tHp = zeros(T,P); %医院t天末报废量
        B_ltp = zeros(35,P,T+1); %t天初交叉配型返回量
        S_tp = zeros(T,P); %t天末短缺量
        vv = zeros(35,P,T); 
        x_tp = zeros(T,P); %调剂量
        y_tpp = zeros(P,P,T); %替代量
        z_tpp = zeros(P,P,T); %调剂替代量
        y_ltp = zeros(35,P,T); %
        C_k = zeros(1,T); %策略成本
        C_EB = zeros(T,P); %血液中心库存成本
        C_EH = zeros(T,P); %医院库存成本
        C_S = zeros(T,P); %短缺成本
        C_Ss = zeros(T,P); %未调剂替代短缺成本
        C_WB = zeros(T,P); %血液中心过期成本
        C_WH = zeros(T,P); %医院过期成本
        C_T = zeros(T,P); %调剂成本
        C_Y = zeros(T,P); %替代成本
        vvv = zeros(35,P,T);
        for t = 1:1:T
            Q_tp(t,:) =  x_kt(t,:)*a_kp; %formula4.1
            for p = 1:1:P
                if (O_tp(t,p) - sum(V_ltHp(1:L_p(p),p,t))) > 0
                    D_tHp(t,p) = (O_tp(t,p) - sum(V_ltHp(1:L_p(p),p,t)));
                else
                    D_tHp(t,p) = 0;
                end
                for l = L_p(p):-1:1
                    if (D_tHp(t,p) - sum(V_ltBp(l:L_p(p),p,t))) > 0
                        q_ltp(l,p,t) = V_ltBp(l,p,t);
                    else
                        if  V_ltBp(l,p,t) - (sum(V_ltBp(l:L_p(p),p,t)) - D_tHp(t,p)) > 0
                            q_ltp(l,p,t) = V_ltBp(l,p,t) - (sum(V_ltBp(l:L_p(p),p,t)) - D_tHp(t,p));
                        else
                            q_ltp(l,p,t) = 0;
                        end
                    end
                end
                for l = L_p(p):-1:1
                    if l == 1
                        V_ltBp(l,p,t+1) = Q_tp(t,p);
                    else
                        V_ltBp(l,p,t+1) = V_ltBp(l-1,p,t) - q_ltp(l-1,p,t);
                    end
                    B_ltp(l,p,t) =  V_ltHp(l,p,t) + q_ltp(l,p,t);
                    if l + gama(p) <= L_p(p)
                        vv(l,p,t) = fix((1-alpha(p))*B_ltp(l,p,t));
                        V_ltHp(l+gama(p),p,t+gama(p)) = vv(l,p,t);
                    end
                end
                if L_p(p) <= 35
                    if V_ltBp(L_p(p),p,t) > q_ltp(L_p(p),p,t)
                        W_tBp(t,p) =  V_ltBp(L_p(p),p,t) - q_ltp(L_p(p),p,t);
                    else
                        W_tBp(t,p) = 0;
                    end
                    if V_ltHp(L_p(p),p,t) > O_tp(t,p)
                        W_tHp(t,p) = V_ltHp(L_p(p),p,t) - O_tp(t,p);
                    else
                        W_tHp(t,p) = 0;
                    end
                end
                if O_tp(t,p) > sum(V_ltHp(:,p,t)) + sum(V_ltBp(:,p,t))
                    S_tp(t,p) = O_tp(t,p) - (sum(V_ltHp(:,p,t)) + sum(V_ltBp(:,p,t)));
                else
                    S_tp(t,p) = 0;
                end
            end
            Ss_tp(t,:) = S_tp(t,:); %未调剂替代短缺量
            for p = P:-1:1
                if S_tp(t,p) > 0 
                    if p ~= 4 && p ~= 5
                        x_tp(t,p) = min(S_tp(t,p),Tt_maxtp(t,p));
                        S_tp(t,p) = S_tp(t,p) - x_tp(t,p);
                        Tt_maxtp(t,p) = Tt_maxtp(t,p) - x_tp(t,p);
                        if S_tp(t,p) > 0
                            for k = 1:1:P
                                if substitue(k,p) ~= 0 && substitue(k,p) ~= inf
                                    if W_tBp(t,k) > 0
                                        y_ltp(L_p(k),k,t) = min(W_tBp(t,k),S_tp(t,p));
                                        W_tBp(t,k) = W_tBp(t,k) - y_ltp(L_p(k),k,t);
                                        S_tp(t,p) = S_tp(t,p) - y_ltp(L_p(k),k,t);
                                        if S_tp(t,p) > 0
                                            for l = L_p(k)-1:-1:1
                                                y_ltp(l,k,t) = min(V_ltBp(l+1,k,t+1),S_tp(t,p));
                                                y_tpp(k,p,t) = y_tpp(k,p,t) + y_ltp(l,k,t) ;                                    
                                                V_ltBp(l+1,k,t+1) = V_ltBp(l+1,k,t+1) - y_ltp(l,k,t);
                                                S_tp(t,p) = S_tp(t,p) - y_ltp(l,k,t);
                                                if S_tp(t,p) == 0
                                                    break
                                                end
                                            end
                                        end
                                    end
                                    if S_tp(t,p) == 0
                                        break
                                    end
                                end
                            end
                            if S_tp(t,p) > 0 
                                for k = 1:1:P
                                    if substitue(k,p) ~= 0 && substitue(k,p) ~= inf
                                        z_tpp(k,p,t) = min(Tt_maxtp(t,k),S_tp(t,p));
                                        Tt_maxtp(t,k) = Tt_maxtp(t,k) - z_tpp(k,p,t);
                                        S_tp(t,p) = S_tp(t,p) - z_tpp(k,p,t);
                                        if S_tp(t,p) == 0
                                            break
                                        end
                                    end
                                end
                            end 
                        end 
                    elseif p == 5
                        if W_tBp(t,4) > 0
                            y_ltp(L_p(4),4,t) = min(W_tBp(t,4),S_tp(t,p));
                            W_tBp(t,4) = W_tBp(t,4) - y_ltp(L_p(4),4,t);
                            S_tp(t,p) = S_tp(t,p) - y_ltp(L_p(4),4,t);
                            if S_tp(t,p) > 0
                                for l = L_p(4)-1:-1:1
                                    y_ltp(l,4,t) = min(V_ltBp(l+1,4,t+1),S_tp(t,p));
                                    y_tpp(4,p,t) = y_tpp(4,p,t) + y_ltp(l,4,t) ;                                    
                                    V_ltBp(l+1,4,t+1) = V_ltBp(l+1,4,t+1) - y_ltp(l,4,t);
                                    S_tp(t,p) = S_tp(t,p) - y_ltp(l,4,t);
                                    if S_tp(t,p) == 0
                                        break
                                    end
                                end 
                            end
                        end
                        if S_tp(t,p) > 0
                            if W_tBp(t,5) > 0
                                y_ltp(L_p(5),5,t) = min(W_tBp(t,5),S_tp(t,p));
                                W_tBp(t,5) = W_tBp(t,5) - y_ltp(L_p(5),5,t);
                                S_tp(t,p) = S_tp(t,p) - y_ltp(L_p(5),5,t);
                                if S_tp(t,p) > 0
                                    for l = L_p(5)-1:-1:1
                                        y_ltp(l,5,t) = min(V_ltBp(l+1,5,t+1),S_tp(t,p));
                                        y_tpp(5,p,t) = y_tpp(5,p,t) + y_ltp(l,5,t) ;                                    
                                        V_ltBp(l+1,5,t+1) = V_ltBp(l+1,5,t+1) - y_ltp(l,5,t);
                                        S_tp(t,p) = S_tp(t,p) - y_ltp(l,5,t);
                                        if S_tp(t,p) == 0
                                            break
                                        end
                                    end 
                                end
                            end
                            if S_tp(t,p) > 0
                                x_tp(t,p) = min(S_tp(t,p),Tt_maxtp(t,p));
                                S_tp(t,p) = S_tp(t,p) - x_tp(t,p);
                                Tt_maxtp(t,p) = Tt_maxtp(t,p) - x_tp(t,p);
                                if S_tp(t,p) > 0
                                    z_tpp(4,p,t) = min(Tt_maxtp(t,4),S_tp(t,p));
                                    S_tp(t,p) = S_tp(t,p) - z_tpp(4,p,t);
                                    Tt_maxtp(t,4) = Tt_maxtp(t,4) - z_tpp(4,p,t);
                                    if S_tp(t,p) > 0
                                        z_tpp(3,p,t) = min(Tt_maxtp(t,3),S_tp(t,p));
                                        S_tp(t,p) = S_tp(t,p) - z_tpp(3,p,t);
                                        Tt_maxtp(t,3) = Tt_maxtp(t,3) - z_tpp(3,p,t);
                                    end
                                end
                            end
                        end         
                    elseif p == 4
                        if W_tBp(t,3) > 0
                            y_ltp(L_p(3),4,t) = min(W_tBp(t,3),S_tp(t,p));
                            W_tBp(t,3) = W_tBp(t,3) - y_ltp(L_p(3),3,t);
                            S_tp(t,p) = S_tp(t,p) - y_ltp(L_p(3),3,t);
                            if S_tp(t,p) > 0
                                for l = L_p(3)-1:-1:1
                                    y_ltp(l,3,t) = min(V_ltBp(l+1,3,t+1),S_tp(t,p));
                                    y_tpp(3,p,t) = y_tpp(3,p,t) + y_ltp(l,3,t) ;                                    
                                    V_ltBp(l+1,3,t+1) = V_ltBp(l+1,3,t+1) - y_ltp(l,3,t);
                                    S_tp(t,p) = S_tp(t,p) - y_ltp(l,3,t);
                                    if S_tp(t,p) == 0
                                        break
                                    end
                                end 
                            end
                        end
                        if S_tp(t,p) > 0
                            x_tp(t,p) = min(S_tp(t,p),Tt_maxtp(t,p));
                            S_tp(t,p) = S_tp(t,p) - x_tp(t,p);
                            Tt_maxtp(t,p) = Tt_maxtp(t,p) - x_tp(t,p);
                            if S_tp(t,p) > 0
                                z_tpp(3,p,t) = min(Tt_maxtp(t,3),S_tp(t,p));
                                Tt_maxtp(t,3) = Tt_maxtp(t,3) - z_tpp(3,p,t);
                                S_tp(t,p) = S_tp(t,p) - z_tpp(3,p,t);
                            end
                        end 
                    end
                end
            end
            C_k(t) = x_kt(t,:)*C_kp; %生成成本
            for p = 1:1:P
                C_EB(t,p) = sum(V_ltBp(:,p,t)) * C_EBp(p);
                C_EH(t,p) = sum(V_ltHp(:,p,t)) * C_EHp(p);
                C_Ss(t,p) = Ss_tp(t,p) * C_Sp(p);
                C_S(t,p) = S_tp(t,p) * C_Sp(p);
                C_WB(t,p) = W_tBp(t,p) * C_Wp(p);
                C_WH(t,p) = W_tHp(t,p) * C_Wp(p);
                C_T(t,p) = (x_tp(t,p) + sum(y_tpp(:,p,t))) *  C_Tp(p);
                for k = 1:1:P
                    if substitue(k,p) ~= 0 && substitue(k,p) ~= inf
                        C_Y(t,p) = (y_tpp(k,p,t) + z_tpp(k,p,t)) * C_ypp(k,p) + C_Y(t,p);
                    end
                end
            end
        end
        for t = 1:1:T
            for p = 1:1:P
                for l = 1:1:L_p(p)
                    VV(l,p,t) =  (V_ltBp(l,p,t) + V_ltHp(l,p,t))*l;
                end
            end
        end
        obj1 = sum(S_tp,'all')/sum(O_tp,'all');   
        obj2 = sum(C_k,'all') + sum(C_EB,'all') + sum(C_EH,'all') +sum(C_S,'all') + sum(C_WB,'all') + sum(C_WH,'all') + sum(C_T,'all') + sum(C_Y,'all');
        obj3 = sum(VV,'all')/(sum(V_ltBp,'all') + sum(V_ltHp,'all'));
        angents_obj(a).obj1 = obj1;
        angents_obj(a).obj2 = obj2;
        angents_obj(a).obj3 = obj3;
    end
end

