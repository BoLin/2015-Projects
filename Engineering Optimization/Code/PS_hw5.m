%Particle Swarm
%x1=[-5 5] x2=[-5 5]
%ploting
x1=linspace(-5,5,100);
x2=linspace(-5,5,100);
[X1,X2]=meshgrid(x1,x2);
f = @(x1,x2) 10*2+(x1.^2-10*cos(2*pi*x1))+(x2.^2-10*cos(2*pi*x2));
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z,0:10:90)
title('Particle Swarm','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
xlabel('X_1','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
ylabel('X_2','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=@(x) 10*2+(x(1)^2-10*cos(2*pi*x(1)))+(x(2)^2-10*cos(2*pi*x(2)));
p=50;%initial particles
X=[0;0];
v(1:2,1:p)=0;%initial v
w=1.4;
c1=2;
c2=2;
dt=1;
loop=1;
Temp=[0;0];

%initial X

for (i=1:p)
    X(1,i)=rand(1)*10-5;
    X(2,i)=rand(1)*10-5;
    v(1,i)=(rand(1)*10-5)/dt;
    v(2,i)=(rand(1)*10-5)/dt;%initial random v 
    i=i+1; 
end
plot(X(1,:),X(2,:),'.')%PLOTING
hold on


Xo=X;
Temp=X(:,1);
for i=1:p
    if f(X(:,i))<f(Temp)
        Temp=X(:,i);
        end
    i=i+1;
end
Xg(1,1:p)=Temp(1,1);
Xg(2,1:p)=Temp(2,1);
%initial iteration and global opt


while loop<200
    r1=rand(1);
    r2=rand(1);

    
    X=X+v*dt;
    
    %updating Xo
        for i=1:p
            if f(Xo(:,i))>f(X(:,i))
             Xo(:,i)=X(:,i);
             end
        i=i+1;
        end
        
        
    
    for i=1:p
        if f(X(:,i))<f(Xg(:,1))
        Xg(1,1:p)=X(1,i);
        Xg(2,1:p)=X(2,i);
        end
    i=i+1;
    end
    
    v=w*v+c1*r1*(Xo-X)/dt+c2*r2*(Xg-X)/dt;
     
    
    for (i=1:p)
   F(i)=f(X(:,i));
        i=i+1;
         end
      
 
    w=1.4*(200-loop)/200;%dynamically 
    loop=loop+1;
    if loop==100
      plot(X(1,:),X(2,:),'+')  
    end
        
end
FX=Xg(:,1)
Fmin=f(Xg(:,1))
plot(X(1,:),X(2,:),'*')%PLOTING
hold on

