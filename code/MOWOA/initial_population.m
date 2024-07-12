function [value] = initial_population(popsize)
    global a_kp T X_1t X_2t
    k = size(a_kp,1);
    i = size(a_kp,2);
    for p = 1:1:popsize
        x_k = zeros(T,k);
        for t = 1:1:T
            for a = 1:1:4
                x_k(t,a) = rand;
            end
            x_k(t,:) = x_k(t,:)./sum(x_k(t,:));
            x_k(t,:) = fix(x_k(t,:)*X_1t(t));
            x_k(t,5) = fix(rand*X_2t(t));
            %x_k(t,6) = X_2t(t) - x_k(t,5);
            x_k(t,6) = fix(rand*X_2t(t));
            if  x_k(t,5) + x_k(t,6) >= X_2t(t)
                x_k(t,6) = X_2t(t) - x_k(t,5);
            end
            Q_tp(t,:) =  x_k(t,:)*a_kp;
            x_k(t,7) = randi([0,Q_tp(t,3)]);
            x_k(t,8) = randi([0,Q_tp(t,3)-x_k(t,7)]);
            x_k(t,9) = randi([0,Q_tp(t,4)+x_k(t,8)]);
            x_k(t,10) = randi([0,Q_tp(t,2)]);
            x_k(t,11) = randi([0,Q_tp(t,6)]);
            Q_tp(t,:) =  x_k(t,:)*a_kp;
            value(p).x_kt = x_k;
        end
    end
end

