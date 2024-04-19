clc;%验证随机序列确实符合拉普拉斯分布
clear;
close;

length=10000;
lambda=0.7;

y=rand_laplace(length,lambda);

[yy,x]=ksdensity(y);

xx=transpose(-5:1e-1:5);
miu=0;
probably=1/(2*lambda)*exp(-abs(xx-miu)/lambda);

figure;
plot(x,yy,'bo');
hold on;
plot(xx,probably,'LineWidth',2);


% plot(xx,probably,'LineWidth',2);
% hold on;
% lambda=2;
% probably=1/(2*lambda)*exp(-abs(xx-miu)/lambda);
% plot(xx,probably,'LineWidth',2);
% xlabel('随机变量X');
% ylabel('概率P');
% legend('b=1','b=2');