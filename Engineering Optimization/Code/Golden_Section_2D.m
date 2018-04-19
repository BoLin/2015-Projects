function[X_min]=Golden_Section_2D(functname,Upperbound,X0)
% optimizing step length parameter using Golden section
% Bracket-> bounds of alpha
% Input X0-> the current design variable vector
% Input S-> search direaction vector
% Output X_min-> minimum of design variables vector
e=0.0001; % Convergence error
N=ceil(log(e)/log(1-0.38197)) % Number of iterations.
Xl=X0;    %Lower bound.
Xu=Upperbound;    %Upper bound.
f=functname;  % objective function.
for i=1:N
    X1=Xl+0.38197*(Xu-Xl);
    X2=Xu-0.38197*(Xu-Xl);
    f1=f(X1);
    f2=f(X2);
    if f1<f2
        Xu=X2;
    elseif f1>f2
        Xl=X1;
    end
    
end
X_min=(Xu+Xl)/2;

        
        
    
    
    