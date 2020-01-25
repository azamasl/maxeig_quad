function [ x,xeigs, X,frec] = SB(pars, options)

% Set options.quad=0 to not use the quadratic approximation
% Set options.quad=1 to starts with line-search method and once dnorm<1e-4
% switch to SBQNKKT method to converge quadratically.
% to obtain dnorm<1e-4 set maxit high enough.
quad = options.quad;
normtol0 = options.normtol;
maxit0 = options.maxit;

options.normtol= 1e-4;
[frec1,y1,Y1,m,dnorm,iter] = SB_ls(pars,options);

options.x0=y1;
options.normtol=normtol0;
options.maxit = maxit0-iter;
fprintf( 'dnorm  = %g .\n',dnorm)

if quad && iter < maxit0 && dnorm < 1e-4
    fprintf( 'USING QUAD approx.\n')
    [ frec2, y2, Y2] = SBQ_NKKT(pars, options,m);
else
    fprintf( '--NOT USING QUAD approx.\n')
    options.maxit = 5; 
    [ frec2, y2, Y2] = SB_ls(pars, options);
end

frec = [frec1(:,1:end-1) frec2];
x= y2;
X = [Y1(:,1:end-1) Y2];
[~,~,~,xeigs]=maxEigSubdiff(x,pars);
