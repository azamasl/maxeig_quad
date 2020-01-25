 function [f, Gnew, m,d,Astry,negC,Q] = maxEigSubdiff(y, pars)
    n=pars.n;
    tol=pars.eigTol;
    Ai= pars.Ai;
    negC=  Ai{1};
    Astry = zeros(size(negC));
    for i=1:n
        Astry = Astry+y(i)*Ai{i+1};
    end

    [Q,D] = eig(Astry+negC);
    d = diag(D);  
    [d,I]=sort(d,'descend');
    Q = Q(:,I);
    lam_1  = d(1);
    
    r = length(d); m=1;  
    if lam_1< 0
        lam_1=0;
        Gnew=zeros(n,1);
    else  
        while m<r && (lam_1-d(m+1))/max(lam_1,1) <= tol
            m=m+1;
        end 
        %This is according to Barvinok-Pataki result.
        m = min(ceil(sqrt(2*n)),m);
        chm = 0.5*(m*(m+1));
        Q1 = Q(:,1:m);
        Gnew = zeros(n,chm);
        for i=1:n
            Gnew(i,:) =mysvec(Q1'*Ai{i+1}*Q1)';
        end

    end
    
 f = lam_1;
end









