function m=fai(h)

e1=1/(2*sqrt(3))*[-2 1 1 1 -2 1];
e2=1/sqrt(6)*[-1 -1 -1 1 1 1];
e3=1/sqrt(6)*[1 1 1 1 1 1];
e4=1/(2*sqrt(7))*[-1 -3 4 -1 1 0];
e5=1/(2*sqrt(21))*[2 -1 -1 -5 -2 7];
e6=1/2*[-1 1 0 -1 1 0];
e=[e1' e2' e3' e4' e5' e6'];
sita=[0 0 0 0 0 0];

for node=1:6
    sita(node)=h'*e(:,node);
end

m=sita;