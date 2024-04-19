clc;
clear;
close all;

%G=graph([1 1 2 2 4 5 3],[2 4 4 3 5 6 6]);
%L=full(laplacian(G));

A=[0.5 0.25 0 0.25 0 0
0.25 0.25 0.25 0.25 0 0
0 0.25 0.5 0 0 0.25
0.25 0.25 0 0.25 0.25 0
0 0 0 0.25 0.5 0.25
0 0 0.25 0 0.25 0.5];

x0=[7 2 3 4 5 6];

n=20000;
state=zeros(n,6);
state(1,:)=x0;
length=size(x0,2);
p=0.15;
noise=[];
aaa=[];

for node=2:n-1
    b=1/(node^p);%生成参数b
    xulie=rand_laplace(length,b);%符合拉普拉斯分布的随机序列
    aaa(node+1,:)=xulie;%储存序列（噪声乘衰减因子等于标准差乘衰减因子）
    m=fai(xulie);%希尔伯特空间的性质，生成随机噪声
    noise(node+1,:)=m;%储存噪声
    state(node,:)=A*(state(node-1,:)+m)';%更新状态
end

e=[];%计算
for node=1:n
    e(node,:)=(state(node,:)-[4.5 4.5 4.5 4.5 4.5 4.5]);
end
    
    
xx=1:1:n;
figure;

plot(xx,state(:,1),'r');%绘制状态图
hold on;
plot(xx,state(:,2),'b');
hold on;
plot(xx,state(:,3),'g');
hold on;
plot(xx,state(:,4),'c');
hold on;
plot(xx,state(:,5),'m');
hold on;
plot(xx,state(:,6),'y');

figure;
plot(xx,noise(:,1));

figure;
plot(xx,e(:,1),'r');%绘制E（x(k)-x的均值）
hold on;
plot(xx,e(:,2),'b');
hold on;
plot(xx,e(:,3),'g');
hold on;
plot(xx,e(:,4),'c');
hold on;
plot(xx,e(:,5),'m');
hold on;
plot(xx,e(:,6),'y');

