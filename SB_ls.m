function [ frec, x, X,m,dnorm,iter] = SB_ls(pars, options)

%Spectral Bundle method with a backtracking line search

x = options.x0;
frec = feval(pars.fgname,x,pars);
subdiffname = pars.subdiffname;
X = x;
dnorm = inf;

n=pars.n;
Ai= pars.Ai;
maxit = options.maxit;

normtol = options.normtol;
gamma = options.gamma;
betta =options.betta;

prtlevel = options.prtlevel;

if prtlevel > 1
    fprintf('\n')
end
for iter = 1:maxit    
   [~, ~,m,dd,Astry,negC,Q] = feval(subdiffname,x, pars);   % maxEigSubdiff %optval = cvx_optval
  
   %% Let the direction to be -grad of the ?Moreau envelope of cutting palne model of the lam_max function
   Q1 = Q(:,1:m);
   Q2 = Q(:,m+1:end);
   P=Q1;
   cvx_begin quiet
        variable Usol(m,m) semidefinite
        W=P*Usol*P';
        expression AW(n,1);
        for i=1:n
            AW(i) = trace(Ai{i+1}'*W);
        end
        maximize trace(P'*(Astry+negC)'*P*Usol) -.5*square_pos(norm(AW))
        subject to 
            trace(Usol)== 1
    cvx_end
    dnew = AW;

    [~,eigsU ] = eig(Usol);
    [dU,~]=sort(diag(eigsU),'descend');
    
    %% Adjusting our guess of the mult. based on the gap in spectrum of Usol
    gm = guessMult(dU, pars.eigTol);
    if(gm ~=m)
        m=gm;    
        Q1 = Q(:,1:m);
        Q2 = Q(:,m+1:end);
        P=Q1;
    end
    
   %% ls  
    t=1;f=feval(pars.fgname,x,pars);   
    fprintf('Iter =%d, t=%g, dnorm: = %g, f = %g, m=%d, eig(1)=%1.16e, eig(2)=%1.16e, eig(3)=%1.16e . \n', iter, t, dnorm, f, m, dd(1), dd(2), dd(3));%5.1e
    dnormnew = norm(dnew);
          
    xnew = x-t*dnew; fnew = feval(pars.fgname,xnew,pars);
    while fnew > f-betta*t*(dnormnew)^2
        t=gamma*t;
        xnew = x-t*dnew; 
        fnew = feval(pars.fgname,xnew,pars);
    end
    
    if dnormnew < normtol
        dnorm = dnormnew;
        x=xnew;
        f=fnew;
        X= [X x];
        frec = [frec f];
        if prtlevel > 0
            fprintf('  tolerance met at iter %d, dnormnew = %g,  f = %g\n', iter,  dnormnew,  f);
        end                
        return 
    end
    if abs(dnormnew-dnorm) < 1e-15
        dnorm = dnormnew;
        x=xnew;
        f=fnew;
        X= [X x];
        frec = [frec f];
        if prtlevel > 0
            fprintf('  dnorm did not decrease sufficiently at iter %d, dnormnew = %g,  dnorm = %g   ,  f = %g\n', iter,  dnormnew,  dnorm, f);
        end                
        return
    end
    dnorm = dnormnew;
    x=xnew;
    f=fnew;
    X= [X x];    
    frec = [frec f];
    
end

end






