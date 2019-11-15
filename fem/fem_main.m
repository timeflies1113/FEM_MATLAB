clear;clc;
%% 零件形状
figure
x_l=[0 0];x_b1=[0 2];x_b3=[4 10];x_r=[10 10];x_t=[10 0];
y_l=[0 -8];y_b1=[-8 -8];y_b3=[-6 -2];y_r=[-2 0];y_t=[0 0];
x_center1=4;y_center1=-8;
n_c=11;%圆弧等分变量n_c 
c1=linspace(pi,pi/2,n_c);r1=2;
x_b2=r1*cos(c1)+x_center1;
y_b2=r1*sin(c1)+y_center1;
x=[x_l,x_b1,x_b2,x_b3,x_r,x_t];
y=[y_l,y_b1,y_b2,y_b3,y_r,y_t];
plot(x,y)
hold on
axis equal
xlabel('x轴');ylabel('y轴'); title('图形'); 


%% 分割区域
tempx=x_t(2);tempy=y_t(2);
x_t(2)=3; y_t(2)=0;%x_b3(1)   x_b1(2)
x_t(3)=1; y_t(3)=0 ;%x_b1(2)    
x_t(4)=tempx; y_t(4)=tempy ;
cut_point=[x_t(2), y_t(2);x_t(3), y_t(3)];    
%分割线数组
x_cutline1=[ x_t(2),x_b3(1)]; 
y_cutline1=[ y_t(2),y_b3(1)];
x_cutline2=[ x_t(3),x_b1(2)]; 
y_cutline2=[ y_t(3),y_b1(2)];
%绘制分割线
plot(x_cutline1, y_cutline1,'m');
plot(x_cutline2, y_cutline2,'m');


%% 边上节点（设置上下边种子节点）
n=[10 30 30 10 10];%各边种子节点数,30*50个单元，31*51个点
x_b1point=linspace(x_b1(1),x_b1(2),n(1)+1);y_b1point=linspace(y_b1(1),y_b1(2),n(1)+1);
x_b3point=linspace(x_b3(1),x_b3(2),n(2)+1);y_b3point=linspace(y_b3(1),y_b3(2),n(2)+1);
x_t1point=linspace(x_t(2),x_t(1),n(3)+1);y_t1point=linspace(y_t(2),y_t(1),n(3)+1);
x_t2point=linspace(x_t(3),x_t(2),n(4)+1);y_t2point=linspace(y_t(3),y_t(2),n(4)+1);
x_t3point=linspace(x_t(4),x_t(3),n(5)+1);y_t3point=linspace(y_t(4),y_t(3),n(5)+1);
%从左至右重排点
x_bpoint=[x_b1point(1:n(1)),x_b2,x_b3point(2:n(2)+1)];y_bpoint=[y_b1point(1:n(1)),y_b2,y_b3point(2:n(2)+1)];
x_tpoint=[x_t3point(1:n(5)),x_t2point,x_t1point(2:n(3)+1)];y_tpoint=[y_t3point(1:n(5)),y_t2point,y_t1point(2:n(3)+1)];


%% 绘制网格、写出各点坐标
hold on
m=30;% 横向网格数目
for i=1:51               
    %画竖线 
    plot([x_tpoint(i),x_bpoint(i)],[y_tpoint(i),y_bpoint(i)],'k');
    % 坐标
    %竖线m等分，直接得到x,y坐标,x1是横坐标,y1是纵坐标，顺序从上之下，从左至右 
    x1(i,:)=linspace(x_tpoint(i),x_bpoint(i),m+1);    
    y1(i,:)=linspace(y_tpoint(i),y_bpoint(i),m+1);
end

% 将坐标排列为一维向量
for j=1:51
    for i=1:31
        x2(i+(j-1)*31)=x1(j,i);
        y2(i+(j-1)*31)=y1(j,i);
    end
end

% 绘制横线
for i=1:m 
    for j=51:-1:2 %j,i与j,i-1点
        plot([x1(j,i) x1(j-1,i)],[y1(j,i) y1(j-1,i)],'k');
    end
end


%% 单元、编号
for j=1:50   %列号，20列单元，21列点
    for i=1:m  %行号，m行单元，m+1行点
        % 从上至下，从左至右编号
        p1=i+(j-1)*(m+1);%第j列，第i行编号（四边形第一个点）
        p2=i+1+(j-1)*(m+1); %第j列，第i+1行编号
        p3=i+1+j*(m+1); %第j+1列，第i+1行编号
        p4=i+j*(m+1);%第j+1列，第i行编号 
        t=i+(j-1)*m; %单元的编号，从上至下，从左至右
        % 单元包括的点的编号
        e(t,:)=[p1,p2,p3,p4];
    end
end


%% 计算总刚度矩阵
E=70E6;%弹性模量
NU=0.25;%泊松比
h=0.02;%厚度
K=zeros(2*31*51,2*31*51);
k=zeros(8,8);
for j=1:30%i是单元的行号，j是单元的列号
    for i=1:50
        t=i+(j-1)*50;
        p1=e(t,1);p2=e(t,2);p3=e(t,3);p4=e(t,4);%提取单元内点的编号到p1,p2,p3,p4
        k=k+stiffness(E,NU,h,x2(p2),y2(p2),x2(p3),y2(p3),x2(p4),y2(p4),x2(p1),y2(p1),1);
        K=assemble(K,k,p2,p3,p4,p1);
    end
end


%% 利用K*u=f求解位移矩阵 改一置零法
% 刚度矩阵K1,将位移为0的对应刚度矩阵对角线置1，其余置0
K1=K;
for i=1:31*2
    K1(i,i)=1;
    for j=i+1:3162
        K1(j,i)=0;
        K1(i,j)=0;
    end
end
% 载荷矩阵f,这是一个集中力载荷(再做一个均布载荷的图形)
f=zeros(3162,1);
f(3102)=-20;
% 解方程，位移矩阵u是每个点u、v的组装
u=K1\f;


%% 变形后的图（位移）
%将u分为x，y方向上的位移
for i=1:31*51
    u1(i)=u(2*i-1);
    u2(i)=u(2*i);
end
% 坐标值加上位移值=位移变形后的图？
for i=1:31*51
    x3(i)=x2(i)+u1(i)*1e7;
    y3(i)=y2(i)+u2(i)*1e7;
end

% 绘制位移图
figure
hold on
% 竖线
for i=1:51
    for j=1:30
        plot([x3(j+(i-1)*31),x3(j+1+(i-1)*31)],[y3(j+(i-1)*31),y3(j+1+(i-1)*31)],'r');
    end
end
% 横线
for i=1:50
    for j=1:31
        plot([x3(j+(i-1)*31),x3(j+i*31)],[y3(j+(i-1)*31),y3(j+i*31)],'r');
    end
end
xlabel('x轴');ylabel('y轴'); title('变形后图形'); 


%% 绘制位移云图(U)
figure
hold on
for i=1:30*50
    j=e(i,:);
    fill(x3(j),y3(j),sqrt(u1(j).^2+u2(j).^2),'FaceColor','interp');
end
shading interp;
colorbar;
colormap jet;
axis equal;
xlabel('x轴');ylabel('y轴'); title('变形后图形上的U位移云图'); 


%% 绘制位移云图(U1)
figure
hold on
for i=1:30*50
    j=e(i,:);
    fill(x3(j),y3(j),u1(j),'FaceColor','interp');
end
shading interp;
colorbar;
colormap jet;
axis equal;
xlabel('x轴');ylabel('y轴'); title('变形后图形上的U1位移云图'); 


%% 绘制位移云图(U2)
figure
hold on
for i=1:30*50
    j=e(i,:);
    fill(x3(j),y3(j),u2(j),'FaceColor','interp');
end
shading interp;
colorbar;
colormap jet;
axis equal;
xlabel('x轴');ylabel('y轴'); title('变形后图形上的U2位移云图'); 



%% RF=K*u;计算反力
RF=K*u;
for i=1:31*51
    RF1(i)=RF(2*i-1);
    RF2(i)=RF(2*i);
end
for i=1:51
    for j=1:31
        RF11(j,i)=RF1(j+(i-1)*31);
        RF22(j,i)=RF2(j+(i-1)*31);
    end
end


%% 绘制反力云图：RF、RF1、RF2
figure
hold on
pcolor(x1',y1',sqrt(RF11.^2+RF22.^2));
colormap jet;
colorbar;
axis equal;
xlabel('x轴');ylabel('y轴'); title('反力RF云图'); 
caxis([0,1]);

figure
hold on
pcolor(x1',y1',RF11);
colormap jet;
colorbar;
axis equal;
xlabel('x轴');ylabel('y轴'); title('反力RF1云图'); 
caxis([0,1]);

figure
hold on
pcolor(x1',y1',RF22);
colormap jet;
colorbar;
axis equal;
xlabel('x轴');ylabel('y轴'); title('支反力RF2云图'); 
caxis([0,1]);


%% 计算单元应力,计算主应力
% 单元应力
stress=zeros(30*50,3);
strain=zeros(30*50,3);
for j=1:30
    for i=1:50
        t=i+(j-1)*50;
        %提取单元内点的编号到p1,p2,p3,p4
        p1=e(t,1);p2=e(t,2);p3=e(t,3);p4=e(t,4);
        %各单元四个点的位移
        uel=[u1(p2);u2(p2);u1(p3);u2(p3);u1(p4);u2(p4);u1(p1);u2(p1)];
        %各单元中心的应力、应变
        [stress(t,:),strain(t,:)]=stressandstrain(E,NU,x2(p2),y2(p2),x2(p3),y2(p3),x2(p4),y2(p4),x2(p1),y2(p1),uel);
    end
end


% 计算主应力大小，方向角
sigma=zeros(30*50,3);
for j=1:30
    for i=1:50
        t=i+(j-1)*50;
        sigma(t,:)=PStresses(stress(t,:));
    end
end


%% 绘制主应力云图前数据处理：节点处应力等于共用该节点单元的应力的平均
com=zeros(31*51,5);%com第一列存储有几个单元共用该点，后四列存储单元编号（不足四列的为0）
for i=1:30*50
    com(e(i,1),1)=com(e(i,1),1)+1;
    com(e(i,1),com(e(i,1),1)+1)=i;
    com(e(i,2),1)=com(e(i,2),1)+1;
    com(e(i,2),com(e(i,2),1)+1)=i;
    com(e(i,3),1)=com(e(i,3),1)+1;
    com(e(i,3),com(e(i,3),1)+1)=i;
    com(e(i,4),1)=com(e(i,4),1)+1;
    com(e(i,4),com(e(i,4),1)+1)=i;
end
%提取sigma第一列
for i=1:50*30
    sigma1(i)=sigma(i,1);
    sigma2(i)=sigma(i,2);
    sigma3(i)=sigma(i,3);
end
% 计算节点处应力值（是周边单元应力的平均值）
sigma1=[sigma1,0];
sigma2=[sigma2,0];
sigma3=[sigma3,0];
for i=1:31*51
    for j=2:5
%         由于调用sigma1时下标不能为0，所以将sigma1后增加一个0，
%         将com(i,j)==0时的sigma1连接它即可解决问题
        if com(i,j)==0 
            com(i,j)=30*50+1;
        end
    end
    %各点的第一主应力
    sigmanode1(i)=(sigma1(com(i,2))+sigma1(com(i,3))+sigma1(com(i,4))+sigma1(com(i,5)))/com(i,1);
    %各点的第二主应力
    sigmanode2(i)=(sigma2(com(i,2))+sigma2(com(i,3))+sigma2(com(i,4))+sigma2(com(i,5)))/com(i,1);
end

%% 绘制第一主应力云图
% 利用fill函数就可以得到四边形的stress云图。利用这个思路，
% 我们在得到有限元计算的所有节点的应力值后，可以对单元进行每
% 个子单元的应力云图绘制，循环完所有的单元后就可以得到整体区域的应力云图了。
figure
hold on
for i=1:30*50
    j=e(i,:);
    fill(x2(j),y2(j),sigmanode1(j),'FaceColor','interp');
end
shading interp;% 去掉网格，使云图更平滑平滑
colorbar;
colormap jet;
axis equal;
xlabel('x轴');ylabel('y轴'); title('第一主应力云图'); 
caxis([-1 20]);


%% 绘制第二主应力云图
figure
hold on
for i=1:30*50
    j=e(i,:);
    fill(x2(j),y2(j),sigmanode2(j),'FaceColor','interp');
end
shading interp;
colorbar;
colormap jet;
axis equal;
xlabel('x轴');ylabel('y轴'); title('第二主应力云图'); 


%% 绘制mise应力、应变云图前的数据处理
for i=1:50*30
    stress1(i)=stress(i,1);
    stress2(i)=stress(i,2);
    stress3(i)=stress(i,3);
    strain1(i)=strain(i,1);
    strain2(i)=strain(i,2);
    strain3(i)=strain(i,3);
end
% 计算节点处应力值（是周边单元应力的平均值）
stress1=[stress1,0];
stress2=[stress2,0];
stress3=[stress3,0];
strain1=[strain1,0];
strain2=[strain2,0];
strain3=[strain3,0];
for i=1:31*51
    for j=2:5
%         由于调用sigma1时下标不能为0，所以将sigma1后增加一个0，
%         将com(i,j)==0时的sigma1连接它即可解决问题
        if com(i,j)==0 
            com(i,j)=30*50+1;
        end
    end
    %各点的stress1平均值
    stressnode1(i)=(stress1(com(i,2))+stress1(com(i,3))+stress1(com(i,4))+stress1(com(i,5)))/com(i,1);
    %各点的stress2平均值
    stressnode2(i)=(stress2(com(i,2))+stress2(com(i,3))+stress2(com(i,4))+stress2(com(i,5)))/com(i,1);
    %各点的stress3平均值
    stressnode3(i)=(stress3(com(i,2))+stress3(com(i,3))+stress3(com(i,4))+stress3(com(i,5)))/com(i,1);
    %各点的mise应力
    mise(i)=(((stressnode1(i)-stressnode2(i))^2+(stressnode2(i)-stressnode3(i))^2+(stressnode3(i)-stressnode1(i))^2)/2)^0.5;
    %各点的strain1平均值
    strainnode1(i)=(strain1(com(i,2))+strain1(com(i,3))+strain1(com(i,4))+strain1(com(i,5)))/com(i,1);
    %各点的strain2平均值
    strainnode2(i)=(strain2(com(i,2))+strain2(com(i,3))+strain2(com(i,4))+strain2(com(i,5)))/com(i,1);
    %各点的strain3平均值
    strainnode3(i)=(strain3(com(i,2))+strain3(com(i,3))+strain3(com(i,4))+strain3(com(i,5)))/com(i,1);
    %各点的应变
    strain0(i)=(((strainnode1(i)-strainnode2(i))^2+(strainnode2(i)-strainnode3(i))^2+(strainnode3(i)-strainnode1(i))^2)/2)^0.5;
end


%% 绘制mise应力云图
figure
hold on
for i=1:30*50
    j=e(i,:);
    fill(x2(j),y2(j),mise(j),'FaceColor','interp');
end
shading interp;% 去掉网格，使云图更平滑平滑
colorbar;
colormap jet;
axis equal;
xlabel('x轴');ylabel('y轴'); title('mise应力云图'); 
caxis([-1 5]);


%% 绘制应变云图
figure
hold on
for i=1:30*50
    j=e(i,:);
    fill(x2(j),y2(j),strain0(j),'FaceColor','interp');
end
shading interp;% 去掉网格，使云图更平滑平滑
colorbar;
colormap jet;
axis equal;
xlabel('x轴');ylabel('y轴'); title('应变云图'); 

% 
% %% pcolor绘图,必须将X,Y,C均设置为网状矩阵
% for i=1:51
%     for j=1:31
%         sigmanode11(j,i)=sigmanode(j+(i-1)*31);
%     end
% end
% figure
% hold on
% pcolor(x1',y1',sigmanode11);
% colormap jet;
% % shading interp;%去掉了网格
% colorbar;
% axis equal;
% xlabel('x轴');ylabel('y轴'); title('sigma1云图'); 
