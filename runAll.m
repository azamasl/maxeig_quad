function runAll(genprob, maxit, eigTol,varargin)

addpath('~/google/projects/mCode/privnlcg');
addpath('~/google/projects/mCode/newcode/subgradmeth');

    n= varargin{1};
    N= varargin{2};       
    if genprob
       pars = genMaxEig(n,N);      
       n=pars.n;
       N=pars.N;
       pars.eigTol = eigTol;
       pars.prob = 'maxEig';
       %% Getting the optimal answer from CVX( uses SDPT3 by defualt)    
       Ai= pars.Ai;       
       tic;
       cvx_begin quiet
          variable Z(N,N) semidefinite
          variable ysol(n,1)
          variable t
          Astary=  Ai{1};
          for i=1:n
            Astary = Astary+ysol(i)*Ai{i+1};
          end
          minimize t
          subject to 
                 Z == t*eye(N) - Astary 
        cvx_end
        CVXtime = toc;
        yopt= ysol;    
        [~,~,~,cvx_eigs]=maxEigSubdiff(ysol, pars);  

    else
       load('data/pars_maxeig', 'pars', 'cvx_eigs','CVXtime', 'yopt');   
    end

    options=setOptions(maxit,pars);       
    %%%x0  = yopt + 1e-4*randn(size(yopt));
    x0= randn(size(yopt));
    options.x0=x0;
    options.quad = 0;

    tic;
    [~,SBlseigs] = SB(pars,options);
    SBlstime=toc;     
    fprintf('--------------------------------------\n');
    options.x0=x0;
    options.quad = 1; 
    tic;
    [~,SBQNKKeigs] = SB(pars,options);        
    SBQNKKtime=toc;       
    Te = table(cvx_eigs, SBlseigs, SBQNKKeigs)  ;
    save(['data/pars_maxeig'], 'options','pars' , 'cvx_eigs','CVXtime', 'yopt');         
end


function options=setOptions( maxit,pars)    
    options.gamma = 0.5;
    options.betta = 1e-4;
    options.evaldist = 1e-9;
    options.normtol = 1e-16;
    options.prtlevel = 0;
    options.maxit = maxit;
    options.x0=abs(randn(pars.nvar,1));   
    options.wolfe1 = 0.0001;
    options.wolfe2 = 0.9;
    options.prtlevel=2;
end  
