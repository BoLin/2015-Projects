%Generalized Reduced Gradient Method
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
%1.3,0.5568,0.3450
X(:,1)=[1.3;0.5568;0.3450];%starting point
alpha=0;%initial alpha value
ii=0;%out iteration
tol=1e-04;


 %functions 
   F =  @(x) x(1)^4-2*x(1)^2*x(2)+x(1)*x(2)^2-2*x(1)+4;
   hh =  @(x) x(1)^2+x(2)^2-2;
   gg =  @(x) 0.25*x(1)^2+0.75*x(2)^2-1+x(3); %adding slack variable x(3)
   

while(1)
    ii=ii+1;%counter
   
    Z=X(1,ii); %Z x1 Y x2 x3
    Y=X(2:3,ii);
   
   
   %calculating gradient
   %finite difference
   GF=gradientf(F,X,ii)'
   GH=[gradientf(hh,X,ii)';gradientf(gg,X,ii)']
   A=[GH(:,1)];
   B=[GH(:,2) GH(:,3)];
   
   GFZ=[GF(:,1)];
   GFY=[GF(:,2);GF(:,3)];
   
   Q=B\A;
   GR=GFZ-Q'*GFY;
   S=-GR
   
   
 %computing stepsize alpha
 q=0;
 
 alpha(1,1)=0;
 x1=X(1,ii);
 x2=X(2,ii);
 
 %x1=(2 - x2^2)^(1/2) h1x
 %x1=(4.0 - 3.0*x2^2)^(1/2) h2x
 
 if abs(x2)>sqrt(2)
     alpha1=-x1/S;
 else
     alpha1=sqrt(2 - x2^2)-x1/S;
 end
 
 if abs(x2)>2*sqrt(3)/3
     alpha2=-x1/S;
 else
     alpha2=sqrt(4.0 - 3.0*x2^2)-x1/S;
 end
 
 alpha(3,1)=min(alpha1,alpha2);
 alpha(2,1)=alpha(3,1)/2;
 
 
 %calculating f(a)
 q=1;%initial q
 while(q<=3)%3 iterations
 DZ=alpha(q)*S;
 Z=Z+DZ;
 DY=-(B\A)*DZ;
    while(1)
    Y=Y+DY;%step 3
    X(:,ii+1)=[Z;Y]
    H=[hh(X(:,ii+1));gg(X(:,ii+1))]
   
      if H<=tol
         break
      else
            DY=B\(-H-A*DZ)
      end
    end%problem with the loop
     fa(q,1)=F(X);%get fa(q)
     q=q+1;
 end
 
 %quadratic interpolation
 x1=alpha(1,1);
 x2=alpha(1,2);
 x3=alpha(1,3);
 f1=fa(1,1);
 f2=fa(2,1);
 f3=fa(3,1);
 
 %computing alphastar
 a2 = (((f3-f1)/(x3-x1))-((f2-f1)/(x2-x1)))/(x3-x2);
 a1 = ((f2-f1)/(x2-x1))-a2*(x1+x2);
 a0 = f1-a1*x1-a2*x1^2;
 FA=@(x) a0+a1*x+a2*x^2;
 alphastar=Golden_Section_2D(FA,x1,x2);
 
 DZ=alphastar*S;
 Z=Z+DZ;
 DY=-(B\A)*DZ;

 while(1)
  Y=Y+DY;%step 3
 X(:,ii)=[Z;Y];
 H=[hh(X(:,ii));gg(X(:,ii))];
    if H==0
       break
    else
        DY=B\-H;
       
    end
 end
 
%Updating X

X(:,ii+1)=[Z;Y]
  
%critieria
DX=X(:,ii+1)-X(:,ii);
if DX'*DX<=tol
    break
end

if ii==50
    break
end

   
end
plot(X(1,:),X(2,:),'LineWidth',2);
