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
    b=1/(node^p);%���ɲ���b
    xulie=rand_laplace(length,b);%����������˹�ֲ����������
    aaa(node+1,:)=xulie;%�������У�������˥�����ӵ��ڱ�׼���˥�����ӣ�
    m=fai(xulie);%ϣ�����ؿռ�����ʣ������������
    noise(node+1,:)=m;%��������
    state(node,:)=A*(state(node-1,:)+m)';%����״̬
end

e=[];%����
for node=1:n
    e(node,:)=(state(node,:)-[4.5 4.5 4.5 4.5 4.5 4.5]);
end
    
    
xx=1:1:n;
figure;

plot(xx,state(:,1),'r');%����״̬ͼ
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
plot(xx,e(:,1),'r');%����E��x(k)-x�ľ�ֵ��
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

