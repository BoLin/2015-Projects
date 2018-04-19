%BFGS method
x1=linspace(-2,2,100);
x2=linspace(-2,2,100);
[X1,X2]=meshgrid(x1,x2);
f =  @(x1,x2)100*((x2-x1.^2)).^2+100*(1-x1).^2;
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z,0:10:500)

hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[0;0]; %starting point

tol=0.0001;
dx=0.01;
a=1.61803;
b=0.0005;
f =  @(x)100*((x(2)-x(1)^2))^2+100*(1-x(1))^2;

n=1;%starting iteration
SIZE=size(x);
H=eye(SIZE(1));%H starting as I matrix

xold=x; %xq-1

while(1)
   x1history(n)=x(1); %save x data to xhistory
   x2history(n)=x(2);
   
   S=@(x) gradientlin(f,dx,xold);
   ST=-H*S(x)' %Search direction
   
   A=0;  %initial A
   minf=@(A) f(x+(A*ST));
   GOLD=GoldSection_1Var(minf,tol,A,a,b,10) 
   A=GOLD(1)
   xold=x;
   x=x+A*ST%update x
   
   p=x-xold;%p is the same
  
  y=gradientlin(f,dx,x)-gradientlin(f,dx,xold);% 3.12b y is y'
  
  sigma=p'*y';%3.22a
  
  t=y*H*y';%3.22b
  
  coe1=(sigma+t)/(sigma.^2);
  coe2=1/sigma;
  UD=coe1*p*p'-coe2*(H*y'*p'+p*(H*y')');%3.20
 
  H=H+UD;%3.19   update H 
 
 E=gradientlin(f,dx,x)-gradientlin(f,dx,xold)
 xold=x;
 
  if norm(E)<=tol
    break
  end
  n=n+1;
  
end
n
history=[x1history' x2history']
plot(x1history,x2history)
axis([-0.5 1.5 -0.5 1.5]); 





