function gradix=gradientlin(functname,dx,x)
format compact
n=size(x);
 N=eye(n(1));
  DeltaX=dx.*N*ones(n(1),1)+[0;0];
  xtemp=x;
for i=1:n;    
    x=xtemp;
    D=[0;0];
    D(i)=DeltaX(i);
    w1=feval(functname,x+D)-feval(functname,x-D);
    w2=2*DeltaX(i);
    gradix(i)=w1/w2;
end

