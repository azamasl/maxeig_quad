function sv = mysvec(S)
    %S is assumed to be symmetric. S_n is isomorphic to R^(n+1;2). 
    %Returns the  column stacking of the lower trianle of S, where 
    %off-diagonal entries are multiplied by sqrt(2)
    s1= size(S,1);
    s2= size(S,2);
    if s1~=s2
       fprintf('S must be square \n');
    end
    C = S-diag(S(eye(s1)==1));
    C = sqrt(2)*C;
    C = diag(S(eye(s1)==1))+ C;    
    %[~,~,sv] = find(mytril(C));
    sv = mytril(C, s2);
end

function sv = mytril(C,s2)
    c=1;
    sv=[];
    for i=1:s2
        sv = [sv;C(c:end,i)];
        c=c+1;
    end
end