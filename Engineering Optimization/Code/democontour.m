x1=linspace(-2,2,100);
x2=linspace(-2,2,100);
[X1,X2]=meshgrid(x1,x2);
f =  @(x1,x2) 3*sin(0.5+2500*x1*x2)*cos(x1);
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z)

x1=linspace(-2,2,100);
x2=linspace(-2,2,100);
[X1,X2]=meshgrid(x1,x2);
f =  @(x1,x2) 3*(sin(0.5+2500*x1*x2))*cos(x1);
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z)
