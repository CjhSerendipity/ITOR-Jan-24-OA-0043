function plotpp_Chinese(pop)
    obj = [pop.obj];
    x = obj(1,:);
    %x = obj(1,:);
    y = obj(2,:);
    z = obj(3,:);
%     [X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'v4');%��ֵ
    %��������
    figure(1)
    %������άɢ��ͼ
    subplot(2,2,1);
    plot3(x,y,z,'*');
    xlabel('��ȱ��');
    ylabel('�ܳɱ�');
    zlabel('ƽ������');
    title('Pareotǰ����');
    grid on;
    %����parotoͼ
%     subplot(2,2,2);
%     surf(X,Y,Z);
    %set(gca,'ztick',0:100000:500000);
    subplot(2,2,2);
    plot(x,y,'*');
    xlabel('��ȱ��');
    ylabel('�ܳɱ�');
    %zlabel('Cost(CNY)');
    %title('Pareot frontier');
    grid on;
    %���������-�ɱ�ͼ
    subplot(2,2,3);
    plot(y,z,'*');
    xlabel('�ܳɱ�');
    ylabel('ƽ������');
    grid on;
    %�������ʶ�-�ɱ�ͼ
    subplot(2,2,4);
    plot(x,z,'*');
    xlabel('��ȱ��');
    ylabel('ƽ������');
    grid on;
 end


