function [f,g,nonsmooth,Astry,negC]=maxEig(y,pars)

%f(y)=\lam_max(A^*y) 
%A^*y = Ai{1} + \sum yi*Ai
%NOTE with C = -Ai{1} we retrieve the standard form of SDP : f(y)=\lam_max(A^*y-C) 
nonsmooth=0;
n=pars.n;
Ai = pars.Ai;
%n=pars.nvar;
negC=  Ai{1};
Astry = zeros(size(negC));
for i=1:n
    Astry = Astry+y(i)*Ai{i+1};
end

[U,D] = eig(Astry+negC);
d= diag(D);
[f,I] = max(d);
if f< 0
    f=0;
    g=zeros(n,1);
else
    U = U(:, I);
    g=[];
    for i=1:n
        g(i) = U'*Ai{i+1}*U;
    end
    g=g';
end
end