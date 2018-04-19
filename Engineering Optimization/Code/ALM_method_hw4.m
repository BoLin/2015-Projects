%Augmented Lagrange Multiplier Method
%ploting
x1=linspace(0,4,100);
x2=linspace(0,4,100);
[X1,X2]=meshgrid(x1,x2);
f = @(x1,x2) x1.^4-2.*(x1.^2).*x2+(x1.^2)+x1.*(x2.^2)-2.*x1+4;
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z,0:1:50)
hold on
ezplot('x1^2+x2^2-2',[0,4,0,4])
hold on
ezplot('0.25*x1^2+0.75*x2^2-1',[0,4,0,4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda1=1;
lambda2=1;
rpmax=1000;
X(:,1)=[3;2]%starting point
rp=rand(1);
gamma=3;
ii=0;%out iteration
while(1)
    ii=ii+1;%counter
    x1=X(1,ii);
    x2=X(2,ii);
   %functions 
   F =  @(x) x(1)^4-2*x(1)^2*x(2)+x(1)*x(2)^2-2*x(1)+4;
   hh =  @(x) x(1)^2+x(2)^2-2;
   gg =  @(x) 0.25*x(1)^2+0.75*x(2)^2-1; 
   AA =  @(x) x(1)^4-2*x(1)^2*x(2)+x(1)*x(2)^2-2*x(1)+4+lambda1*( x(1)^2+x(2)^2-2)+lambda2*( 0.25*x(1)^2+0.75*x(2)^2-1);
  %minimize A(x,lambda,rp)
  X=BFGSfunction(AA,X,ii);
  
  %convergence criterion    
 if rp>rpmax
    break
 end
 
 ff(ii)=F(X(:,ii));
 if ii>1
    if ff(ii)-ff(ii-1)<1e-04
    break
    end
end
%f(i)=h(X(:,i));
 gv(ii)=gg(X(:,ii));
 hv(ii)=hh(X(:,ii));
    
%update iteration    
    lambda1=lambda1+2*rp*max(gv(ii),-lambda1/(2*rp));
    lambda2=lambda2+2*rp*hv(ii);
    rp=rp*gamma;
end
plot(X(1,:),X(2,:),'LineWidth',2);


