function x=rand_laplace(siz,lambda)

% if nargin==1%�����������Ŀ
%     lambda=sqrt(2);%���û��lambda�͸�һ��ֵ
% end

z=rand(siz,1);
x=zeros(siz,1);
% lambda=2;
in=z<=0.5;
ip=z>0.5;
x(in)=lambda*log(2*z(in));
x(ip)=-lambda*log(2*(1-z(ip)));
