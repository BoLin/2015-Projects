%Midtern Genetic algorithm
%x1=[0 5] x2=[0 8]
%ploting
x1=linspace(0,5,100);
x2=linspace(0,8,100);
[X1,X2]=meshgrid(x1,x2);
f = @(x1,x2) 3*sin(0.5+0.25*x1.*x2).*cos(x1);
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z)
title('Genetic Algorithm','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
xlabel('X_1','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
ylabel('X_2','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Midtern Method of Feasible Direction
%x1=[0 3] x2=[0 2]
%ploting
x1=linspace(0,3,100);
x2=linspace(0,2,100);
[X1,X2]=meshgrid(x1,x2);
f = @(x1,x2) -x1.*x2;
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z)
hold on
ezplot('(x1^2)/9+(x2^2)/4-1',[0,3,0,2])
title('Method of Feasible Direction','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
xlabel('X_1','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
ylabel('X_2','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = @(x1,x2) -x1.*x2;
g = @(x1,x2) (x1.^2)/9+(x2.^2)/4-1 
x1=2.7529
x2=0.75826
df_dx1 = (f(1.001*x1,x2)-f(x1,x2))/(0.001*x1);
df_dx2 = (f(x1,1.001*x2)-f(x1,x2))/(0.001*x2);
df=[df_dx1;df_dx2]
S=-df
X=[0;0]

%[a]=solve('((x1+a*S(1))^2)/9+((x2+a*S(2))^2)/4-1=0')
[a]=solve('((2.7529-a*0.91617)^2)/9+((0.75826+a*0.4008)^2)/4-1=0')
X(:,1)=[1 1]
X(:,2)=[1.6641 1.6641]
X(:,3)=[2.7529 0.75826]
X(:,4)=[-0.0838 1.9992]
T= @(x1,x2) -x1*x2+1*10*((((x1^2)/9+(x2^2)/4-1)/-0.1)^2-3*((x1^2)/9+(x2^2)/4-1)+3)

plot(X(1,:),X(2,:),'LineWidth',2)

dg_dx1 = (g(1.001*x1,x2)-g(x1,x2))/(0.001*x1);
dg_dx2 = (g(x1,1.001*x2)-g(x1,x2))/(0.001*x2);
dg=[dg_dx1;dg_dx2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=linspace(0,3,100);
x2=linspace(0,2,100);
[X1,X2]=meshgrid(x1,x2);
T= @(x1,x2) -x1.*x2+1*1*10*((((x1.^2)/9+(x2.^2)/4-1)/-0.1).^2-3*((x1.^2)/9+(x2.^2)/4-1)+3)
f(X1,X2);
Z=T(X1,X2);
contourf(X1,X2,Z,30)
hold on
ezplot('(x1^2)/9+(x2^2)/4-1',[0,3,0,2])
title('r_p=1','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
xlabel('X_1','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
ylabel('X_2','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
hold on