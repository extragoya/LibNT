function is_eq = fuzzy_eq(A,B,tol)
%EQ Summary of this function goes here
%   Detailed explanation goes here

    
    function is_eq=fuzzy_compare(a,b,tol)
       test=abs(a-b);
       inds=find(test>tol, 1);
       if(isempty(inds))
           is_eq=true;
       else
           is_eq=false;
       end
    end
    func=@(a,b)fuzzy_compare(a,b,tol);
    is_eq= compare(A,B,func);
end

