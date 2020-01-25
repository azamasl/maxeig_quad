function [ frec, x, X] = SBQ_NKKT(pars, options,mopt)

% Using "Null-Space" method to solve the KKT condition arising from applying
% newton methood to max_eig problem (See page 531 of N&W book). 
% For more info refer to the following papers
%[HO]: "The spectral bundle method with second-order information",
% C. Helmberg, M. L. Overton, F. Rendl
%[OW]: "Second Derivatives for Optimizing Eigenvalues of Symmetric
%Matrices", Michael L. Overton, Robert S. Womersley.

% x0 must be in a close neigbourhood of x*, 

fprintf('SBQ_NKKT.\n')

maxit = options.maxit;
normtol = options.normtol;
prtlevel = options.prtlevel;
x = options.x0;
X = x;
frec = feval(pars.fgname,x,pars);
subdiffname = pars.subdiffname;

%primal-dual variable:
pdvar = [abs(randn(1,1)); x];
dnorm = inf;
n=pars.n;
N=pars.N;
Ai= pars.Ai;
AbarT = [];%N^2 by n matrix
for i=1:n
    AbarT = [AbarT vec(Ai{i+1})];
end

HUtilda =eye(n);

for iter = 1:maxit    

    [f, ~,~,dd,~,~,Q] = feval(subdiffname, pdvar(2:end), pars); 
    m=mopt;    chm = 0.5*(m*(m+1));
    Q1 = Q(:,1:m);  
    Q2 = Q(:,m+1:end);
    %Util2 is the minimizing subgrad
    cvx_begin  quiet
        cvx_precision high
        variable Util2(m,m) symmetric
        W=Q1*Util2*Q1';
        expression AW(n,1);
        for i=1:n
            AW(i) = trace(Ai{i+1}*W);
        end
        minimize dot(AW,AW)
        subject to 
            trace(Util2)== 1
    cvx_end
    
    % Compute HUtilda (Hessian) from  (9) in [HO].   
    lam1 = max(dd);
    D2 = diag(lam1-dd(m+1:end));%of size N-m
    D2inv = inv(D2);
    HUtil2 = zeros(n,n);
    for i=1:n
        for j=1:n
            HUtil2(i,j) = 2*trace(Ai{i+1}*Q1*Util2*Q1'*Ai{j+1}*Q2*D2inv*Q2');
        end
    end  

    % Construct b_hat K_hat from (4.5),(4.6)&(4.7) in [OW]
    b_hat = mysvec(diag(dd(1:m)-lam1));%chm x 1    
    K_hat = zeros(chm,n+1);
    K_hat(:,1)  = mysvec(eye(m));
    for i=2:n+1
        K_hat(:,i) =-mysvec(Q1'*Ai{i}*Q1);
    end
    
    K_hatT = K_hat';% n+1 by chm 
    %QK:n+1 x n+1
    %RK:%n+1 x chm    
    [QK, RK] = qr(K_hatT);
    Y = QK(:,1:chm);%Y is the range of K_hat'
    Z = QK(:,chm+1:end);
    %assuming chm<n+1,     
    RK = RK(1:chm,:) ;%chm x chm
    p = RK'\b_hat;%chm x 1
    Yp = Y*p;
    
    W = [0  zeros(1,n); zeros(n,1)  HUtil2];
    ZtWZ = Z'*W*Z;
    bb = -Z'*([1;zeros(n,1)] + W*Yp);
    q = ZtWZ\bb;
    deltapdvar = QK*[p;q];     
    pdvar = pdvar + deltapdvar;
    
    %% updating x and f
    fprintf('--------Iter =%d,  dnorm: = %g, f = %g, eig(1)=%1.16e, eig(2)=%1.16e, eig(3)=%1.16e . \n', iter, dnorm, f,  dd(1), dd(2), dd(3));%5.1e 
    x = pdvar(2:end);    
    f = feval(pars.fgname,x,pars);           
    X= [X x];    
    frec = [frec f]; 
        
    dnormnew = norm( deltapdvar(2:end));
    
    if dnormnew < normtol   
        dnorm = dnormnew; 
        if prtlevel > 0
            fprintf('  tolerance met at iter %d, dnormnew = %g,  f = %g\n', iter,  dnormnew,  f);
        end                
        return 
    end
    if abs(dnormnew-dnorm) < 1e-15
        dnorm = dnormnew; 
        if prtlevel > 0
            fprintf('  dnorm did not decrease sufficiently at iter %d, dnormnew = %g,  dnorm = %g   ,  f = %g\n', iter,  dnormnew,  dnorm, f);
        end                
        return
   end 
   dnorm = dnormnew;  
end

end




