function c = guessMult(d, tol)
   lam1 = max(d);
   c=0;
   for i = 1:length(d)
       if abs(d(i)-lam1)/max(lam1,1) <= tol
            c=c+1; 
       end
   end 
   
end