%SQP Method
%ploting
% x1=linspace(0,4,100);
% x2=linspace(0,4,100);
% [X1,X2]=meshgrid(x1,x2);
% f = @(x1,x2) x1.^4-2.*(x1.^2).*x2+(x1.^2)+x1.*(x2.^2)-2.*x1+4;
% f(X1,X2);
% Z=f(X1,X2);
% contour(X1,X2,Z,0:1:50)
% hold on
% ezplot('x1^2+x2^2-2',[0,4,0,4])
% hold on
% ezplot('0.25*x1^2+0.75*x2^2-1',[0,4,0,4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X(:,1)=[3;2;0.3450];%starting point
alpha=0;%initial alpha value
ii=0;%out iteration
tol=1e-04;
H=eye(2)%initial hessian
alpha=0;%initial alpha
   
F =  @(x) x(1)^4-2*x(1)^2*x(2)+x(1)*x(2)^2-2*x(1)+4;
hh =  @(x) x(1)^2+x(2)^2-2;
gg =  @(x) 0.25*x(1)^2+0.75*x(2)^2-1;

while(1)
 ii=ii+1;   

 S=quadprog(H,F,A,b,Aeq,beq);%QP function output S
%calculate lambda

%calculating alpha
phy=@(x) x(1)^4-2*x(1)^2*x(2)+x(1)*x(2)^2-2*x(1)+4+lambda(j)*(x(1)^2+x(2)^2-2)+lambda(m+k)*(max((0.25*x(1)^2+0.75*x(2)^2-1),0));
minf=@(A) phy(x+(alpha*S));
alpha=GoldSection_1Var(minf,tol,A,a,b,10) 


DX=alpha*S;
X(:,ii)=X(:,ii-1)+DX;
%convergence for SQP
 if hh(X)<tol && gg(X)<tol
    break
      
 end
 
 
 
 
 %critieria

if DX'*DX<=tol
    break
end

if ii==50
    break
end
%update Metric H
y=DX;
p=
H=H+(p*p'/p'*y)-(H'*H)/(y'*H*y);
 
 
 
 end
 

