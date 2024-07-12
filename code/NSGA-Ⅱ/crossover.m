function [newChromo1,newChromo2] = crossover(chromo1,chromo2,pc)
    global T;
    newChromo1 = chromo1;
    newChromo2 = chromo2;
    n = rand;
    if n<=pc
        for t = 1:1:T
            a = randi(11);
            b = randi(11);
            newChromo1((t-1)*11+1:(t-1)*11+a) = chromo2((t-1)*11+1:(t-1)*11+a);
            newChromo2((t-1)*11+1:(t-1)*11+b) = chromo1((t-1)*11+1:(t-1)*11+b);
        end
    end
%         %第一段交叉
%         for i=1:1:iNum
%             a = randi(jNum);%任选一点进行交叉
%             b = randi(jNum);
%             newChromo1S1((i-1)*jNum+1:(i-1)*jNum+a) = chromo2S1((i-1)*jNum+1:(i-1)*jNum+a);
%             newChromo2S1((i-1)*jNum+1:(i-1)*jNum+a) = chromo1S1((i-1)*jNum+1:(i-1)*jNum+a);
%             newChromo1S2((i-1)*jNum+1:(i-1)*jNum+b) = chromo2S2((i-1)*jNum+1:(i-1)*jNum+b);
%             newChromo2S2((i-1)*jNum+1:(i-1)*jNum+b) = chromo1S2((i-1)*jNum+1:(i-1)*jNum+b);
%         end
%         %第二段交叉
%         for i=1:1:iNum
%             newChromo1S1(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i) = fix(chromo1S1(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.7+chromo2S1(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.3);
%             newChromo2S1(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i) = fix(chromo1S1(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.3+chromo2S1(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.7);
%             newChromo1S2(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i) = fix(chromo1S2(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.7+chromo2S2(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.3);
%             newChromo2S2(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i) = fix(chromo1S2(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.3+chromo2S2(iNum*jNum+(i-1)*jNum+1:iNum*jNum+jNum*i).*0.7);
%         end
%     end
%     newChromo1 = [newChromo1S1,newChromo1S2];
%     newChromo2 = [newChromo2S1,newChromo2S2];
end
