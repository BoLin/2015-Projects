function[Bracket]=Bracketing_2D_HW3(X0,S)
% finding step length bracket or one-dimensional parameter bracket
% Input X0-> the current design variable vector
% Input S-> search direaction vector
% Output Bracket-> lower and upper bounds for alpha
al=0;
au_max=min(([5;5]-abs(X0))./abs(S));  % Uni-modal upper bound.

f=@(x1,x2) 100*(x2-x1^2).^2+100*(1-x1)^2;  % objective function.
au=al+0.01*(au_max-al);
i=0;
while au < au_max
    
        a1=au;  
        au=(2.6154)*a1-1.6154*al;
        X1=X0+a1*S;
        Xu=X0+au*S;
        if f(Xu(1),Xu(2))>f(X1(1),X1(2))
                Bracket=[al,min(au,au_max)];
                break;
        end
        
         al=a1;
end
