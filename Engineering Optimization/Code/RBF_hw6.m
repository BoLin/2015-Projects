%HW6_1
%Radial basis function


%LHS 15 samples
x = lhsdesign(16,2)%16 samples 2 varaibles
X=[x(:,1)*5 x(:,2)*8]
tol=0.1


for i=1:16
    Y(i,1)= 3*sin(0.5+0.25*X(i,1)*X(i,2))*cos(X(i,1))
end

% create a neural network
net = newrb(X',Y')

%initial 16 parents
p=16;%number of parents
ns=16;%number of designs
c=5;%not -3
Sumf=0;%initial sum of fitness
loop=0;
Frank=[0;0];

%initial X
X=X';

f = @(x) 3*sin(0.5+0.25*x(1)*x(2))*cos(x(1));
while (1)
% create a neural network
for i=1:16
    Y(i,1)= 3*sin(0.5+0.25*X(1,i)*X(2,i))*cos(X(1,i))
end
net = newrb(X,Y')
    
    
   Sumf=0; 
      for (i=1:p)
   F(i)=sim(net,[X(1,i);X(2,i)]);%calculating n(x)
    Fs(i)=c-F(i);%minimization fitness
    Sumf=Sumf+Fs(i);
    i=i+1;
         end
    

F%value of fx

%rank based selection
    
% end of GA situation
    for (i=1:p)
    rank=p;
            for (j=1:p)
                 if Fs(i)<Fs(j)
                 rank=rank-1;
                 end
            j=j+1;
            end
   Frank(1,i)=rank;
   Frank(2,i)=Fs(i);%Frank is with X
   i=i+1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%
    



if loop<=80
    r=Fs/Sumf;%5.4.7
else

    for(i=1:p)
r(i)=2*(ns-i+1)/((ns+1)*ns);%5.4.11
    end
end
    
    Ran=rand(1,p);%p random fall


Rf(1)=0;%wheel initialize
for (i=2:p+1)
    Rf(i)=Rf(i-1)+r(i-1);
    i=i+1;
end
%normal parrent choosing
for (i=1:p)
    ii=1;
    while ii<=p
        if Ran(i)<Rf(ii+1)
            Ran1(i)=ii;
            ii=p+15;
        else
        ii=ii+1;
        end
    end
    i=i+1;
end


%special transfer for ranking fitness
%if loop>80%same
for (i=1:p)
    for (j=1:p)
     if Ran1(i)==Frank(1,j)
         Ran1(i)=j;%Ran1 is erased
     end
     j=j+1;
    end
    i=i+1;
end
%end
X=(X+5)*51.1%DEC to BIN
%binary lenth 9 (x+5)*51.1=new x for dec2bin
k=randi([1,8]);%random string length
%mutation
for i=1:p-1
 sq1=Ran1(i);
 sq2=Ran1(i+1);
[C(:,i),C(:,i+1)]=Mute(X(:,sq1),X(:,sq2),k);
i=i+2;
end

%Elitist Strategy

T=[X C];
T=(T/51.1)-5;%back to dec
for i=1:2*p
    for j=i:2*p
        if sim(net,[T(1,i);T(2,i)])>sim(net,[T(1,j);T(2,j)])
            TV=T(:,j);
            T(:,j)=T(:,i);
            T(:,i)=TV;
        end
        j=j+1;
    end
    i=i+1;
end

X=T(:,1:p);


plot(X(1,:),X(2,:),'o')%PLOTING
hold on

if  abs(f(X(:,1))+3)<tol
    break
end

loop=loop+1;
end

 Fmin=f(X(:,1))%the min f
 X(:,1)
 plot(X(1,2),X(2,2),'*')
 hold on
 loop

 
 
x1=linspace(0,5,16);
x2=linspace(0,8,16);
[XF1,XF2]= meshgrid(x1,x2)


for i=1:16
    for j=1:16
        N(i,j)=sim(net,[XF1(i,j);XF2(i,j)]);
    end
end

contour(XF1,XF2,N)
title('RDF','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
xlabel('X_1','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
ylabel('X_2','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
hold on
axis([0 5 0 8])




% 
% x1=linspace(0,5,100);
% x2=linspace(0,8,100);
% [X1,X2]=meshgrid(x1,x2);
% f = @(x1,x2) 3*sin(0.5+0.25*x1.*x2).*cos(x1);
% 
% f(X1,X2);
% Z=f(X1,X2);
% contour(X1,X2,Z)
% title('The Original Model','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
% xlabel('X_1','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
% ylabel('X_2','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
% hold on