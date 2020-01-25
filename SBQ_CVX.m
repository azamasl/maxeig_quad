function [ frec, x, X] = SBQ_CVX(pars, options,mopt)

% It uses CVX to get both the direction as well as
% to approximate the Hessian. See SBQ_NKKT for a comparison, which only
% uses CVX to compute the Hessian. It then solves the linear system exactly. 
% x0 must be in a close neigbourhood of x*, 

fprintf('SBQ_CVX. \n')

maxit = options.maxit;
normtol = options.normtol;
prtlevel = options.prtlevel;
x = options.x0;
X = x;
frec = feval(pars.fgname,x,pars);
subdiffname = pars.subdiffname;

dnorm = inf;
n=pars.n;
N=pars.N;
Ai= pars.Ai;
AbarT = [];%needs to be N^2 by n matrix
for i=1:n
    AbarT = [AbarT vec(Ai{i+1})];
end


HUtilda =eye(n);

for iter = 1:maxit    
    
   % Let direction to be neg grad of the moreau env of cutting palne modle of lam_max fun
   [f, ~,~,dd,Astry,negC,Q] = feval(subdiffname,x, pars);   % maxEigSubdiff %optval = cvx_optval   
   m=mopt;
   Q1 = Q(:,1:m);
   Q2 = Q(:,m+1:end);
   P=Q1;   
   
   %% Estimating the Hessian to get the dir.    
   cvx_begin  quiet
        cvx_precision high
        variable Utilda(m,m) symmetric
        W=P*Utilda*P';
        expression AW(n,1);
        for i=1:n
            AW(i) = trace(Ai{i+1}'*W);
        end
        minimize dot(AW,AW)
        subject to 
            trace(Utilda)== 1
    cvx_end
    [~,pos] = chol(Utilda);
    if pos
        fprintf('Utilda is not SPD \n');
    end 
    
    %I'm being careful here:
    Util = zeros(size(Utilda));
    for i = 1:size(Utilda,1)
        for j = 1:size(Utilda,2)
            Util(i,j) = Utilda(i,j);
        end
    end
    
    Qkron = kron( Q1,Q2);
    lam1 = max(dd);
    D2 = diag(lam1-dd(m+1:end));%should be of size N-m
    D2inv = inv(D2);
    HUtilda = kron(Util,D2inv);
    HUtilda =   Qkron*HUtilda*Qkron';
    HUtilda = 2*AbarT'*HUtilda*AbarT; 

    %% Using CVX to get the direction % NOT relaxing the equality const. in OW4
    cvx_begin quiet
        cvx_precision best
        variable yq(n,1) 
        variable delta
        Astaryq= negC;
        for i=1:n
            Astaryq = Astaryq+yq(i)*Ai{i+1};
        end
        minimize trace((yq-x)'*HUtilda'*(yq-x))+ delta        
        subject to 
            P'*Astaryq*P== delta*eye(m)
   cvx_end
   dnew = yq-x;   
   %% updating x and f 
   x=yq;
   f = feval(pars.fgname,x,pars);
   X= [X x];    
   frec = [frec f]; 
   
   dnormnew = norm(dnew);
   fprintf('Iter =%d,  dnorm: = %g, f = %g, eig(1)=%1.16e, eig(2)=%1.16e, eig(3)=%1.16e . \n', iter, dnormnew, f,  dd(1), dd(2), dd(3));%5.1e       
   if dnormnew < normtol   
        %frec = [frec f];
        if prtlevel > 0
            fprintf('  tolerance met at iter %d, dnormnew = %g,  f = %g\n', iter,  dnormnew,  f);
        end                
        return 
   end
   if abs(dnormnew-dnorm) < 1e-15
        if prtlevel > 0
            fprintf('  dnorm did not decrease sufficiently at iter %d, dnormnew = %g,  dnorm = %g   ,  f = %g\n', iter,  dnormnew,  dnorm, f);
        end                
        return
   end 
   dnorm = dnormnew;
end

end