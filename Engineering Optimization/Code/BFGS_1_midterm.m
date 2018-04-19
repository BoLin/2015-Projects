%BFGS MIDTERM
f = @(x1,x2)3*(sin(0.5+2500*x1*x2))*cos(x1);
ezcontour(f,[0,5],[0,0.001])%axis([0 5 -0.001 0.005]); 
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=[0.5;0.0007]; %starting point
tol=0.0001;
dx=0.01;
a=1.61803;
b=0.0005;
f = @(x) 3*(sin(0.5+2500*x(1)*x(2)))*cos(x(1));

n=1;%starting iteration
SIZE=size(x);
H=eye(SIZE(1));%H starting as I matrix

xold=x; %xq-1
ni=0;

while(n<200)
    nr=0
   x1history(n)=x(1); %save x data to xhistory
   x2history(n)=x(2);
   
   S=@(x) gradientlin(f,dx,xold);
   ST=-H*S(x)' %Search direction
   
   A=0;  %initial A
   minf=@(A) f(x+(A*ST));
   GOLD=GoldSection_1Var(minf,tol,A,a,b,10) 
   A=GOLD(1)
   xold=x
   xold2=x;
   x=x+A*ST%update x
   ftemp=feval(f,x)
   
   
    while x(2)<0 || x(2)>0.0008||x(1)<0.5||x(1)>5
                fprintf('out of bound iteration\n')
        x=xold
        HR=[sqrt(3)/2 -1/2;1/2 sqrt(3)/2];%rotate search direction
        %HR=[0 -1;1 0];
        x=x+A*HR*ST
        ST=HR*ST;
        if nr==11
         A=0.9*A;
         % A=0.9*A;
          nr=0;
        end 
       nr=nr+1;
    end
      
   
   p=x-xold;%p is the same
  
  y=gradientlin(f,dx,x)-gradientlin(f,dx,xold);% 3.12b y is y'
  
  sigma=p'*y';%3.22a
  
  t=y*H*y';%3.22b
  
  coe1=(sigma+t)/(sigma.^2);
  coe2=1/sigma;
  UD=coe1*p*p'-coe2*(H*y'*p'+p*(H*y')');%3.20
 H
 
 if ni==3
   H=H+UD%3.19   update H 
  ni=0;
 end
  
 
 %E=gradientlin(f,dx,x)-gradientlin(f,dx,xold)
 E=gradientlin(f,dx,x)
 
 
  if norm(E)<=1.0e-004
     condition=1;
      break
  end
  
  if abs((feval(f,x)-feval(f,xold)))<=0.0001
      condition=2;
      %break
  end
   
  
  n=n+1;
  ni=ni+1;
 
end
n
condition
history=[x1history' x2history']
plot(x1history,x2history,'LineWidth',1)
axis([0 5 0 0.001]); 