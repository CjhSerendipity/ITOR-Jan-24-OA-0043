function plotpp_Chinese(pop)
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
    xlabel('短缺率');
    ylabel('总成本');
    zlabel('平均库龄');
    title('Pareot前沿面');
    grid on;
    %绘制paroto图
%     subplot(2,2,2);
%     surf(X,Y,Z);
    %set(gca,'ztick',0:100000:500000);
    subplot(2,2,2);
    plot(x,y,'*');
    xlabel('短缺率');
    ylabel('总成本');
    %zlabel('Cost(CNY)');
    %title('Pareot frontier');
    grid on;
    %绘制满意度-成本图
    subplot(2,2,3);
    plot(y,z,'*');
    xlabel('总成本');
    ylabel('平均库龄');
    grid on;
    %绘制新鲜度-成本图
    subplot(2,2,4);
    plot(x,z,'*');
    xlabel('短缺率');
    ylabel('平均库龄');
    grid on;
 end


