x1=linspace(-2,2,100);
x2=linspace(-2,2,100);
[X1,X2]=meshgrid(x1,x2);
f =  @(x1,x2)100*((x2-x1.^2)).^2+100*(1-x1).^2;
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z,0:10:500)

hold on

x=[0;0];


tol=0.0001;
dx=0.01;
a=1.61803;
b=0.0005;
f =  @(x)100*((x(2)-x(1)^2))^2+100*(1-x(1))^2;

n=1;

while(1)
  x1history(n)=x(1);
  x2history(n)=x(2);
  A=0;
 S=@(x) -1*gradientlin(f,dx,x);
 ST=S(x)'
 minf=@(A) f(x+(A*ST));
 GOLD=GoldSection_1Var(minf,tol,A,a,b,10)
 
 A=GOLD(1)
 %feval(f,x)
 %feval(f,x+(A*ST))
 E=feval(f,x)-feval(f,x+(A*ST))
 p=A*ST
 x=x+p
  if abs(E)<=tol
    break
  end
  n=n+1;
  
end
n
history=[x1history' x2history']
plot(x1history,x2history)
axis([-0.5 1.5 -0.5 1.5]); 





