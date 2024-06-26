clc;
clear;
close all;

%G=graph([1 1 2 2 4 5 3],[2 4 4 3 5 6 6]);
%L=full(laplacian(G));

A=0.25*[2 1 0 0 1
1 2 1 0 0
0 1 2 0 1
0 0 0 3 1
1 0 1 1 1];

x0=[7 2 3 4 5];

n=50;
state=zeros(n,5);
state(1,:)=x0;
length=size(x0,2);
p=0.15;
noise=[];
aaa=[];

for node=1:n%生成正态分布序列
    x=normal(5);
    aaa(node,:)=x;
end

fa=0.9;

for node=2:n%生成随机噪声
    noise(node+1,:)=fa^node*aaa(node,:)-fa^(node-1)*aaa(node-1,:);
end

for node=2:n-1%更新状态
    state(node,:)=A*(state(node-1,:)+noise(node-1,:))';
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

figure;
plot(xx,noise(1:50,1));
