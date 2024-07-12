function plotpp(pop)
    obj = [pop.obj];
    x = obj(1,:);
    %x = obj(1,:);
    y = obj(2,:);
    z = obj(3,:);
%     [X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'v4');%插值
    %创建窗口
    figure(1)
    %绘制三维散点图
    subplot(2,2,1);
    plot3(x,y,z,'*');
    xlabel('Shortage rate');
    ylabel('Cost(CNY)');
    zlabel('Average stock age');
    title('Pareot frontier');
    grid on;
    %绘制paroto图
%     subplot(2,2,2);
%     surf(X,Y,Z);
    %set(gca,'ztick',0:100000:500000);
    subplot(2,2,2);
    plot(x,y,'*');
    xlabel('Shortage rate');
    ylabel('Cost(CNY)');
    %zlabel('Cost(CNY)');
    %title('Pareot frontier');
    grid on;
    %绘制满意度-成本图
    subplot(2,2,3);
    plot(y,z,'*');
    xlabel('Cost(CNY)');
    ylabel('Average stock age');
    grid on;
    %绘制新鲜度-成本图
    subplot(2,2,4);
    plot(x,z,'*');
    xlabel('Shortage rate');
    ylabel('Average stock age');
    grid on;
 end


