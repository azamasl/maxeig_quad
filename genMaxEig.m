function pars = genMaxEig(n,N)     
    fprintf('Creating a random maxEig problem \n');
   
    Ai={};
    Brec=[];
    for i=1:n+1        
        B =randn(N,N);        
        Ai{i}=(B+B')./2;
        Brec = [Brec ; B];
    end
   
    fid = fopen('/Users/azam/google/projects/eclipsews/CBundle/eig.txt','w');
    fprintf(fid, '%g ', n);
    fprintf(fid, '%g\n', N);
    l=N*(n+1);
    for i=1:l
        fprintf(fid, '%1.15f ', Brec(i,1:end-1));
        fprintf(fid, '%1.15f', Brec(i,end));
        fprintf(fid, '\n');
    end
    fprintf(fid, '%1.15f ', Brec(n+1,1:end-1));
    fprintf(fid, '%1.15f', Brec(n+1,end));
    fclose(fid);
    
    
    pars.n=n;
    pars.N= N;
    pars.Ai = Ai;
    pars.nvar =n; 
    pars.fgname = 'maxEig';
    pars.subdiffname = 'maxEigSubdiff';
end
