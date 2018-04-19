%Genetic algorithm
%x1=[-5 5] x2=[-5 5]
%ploting
x1=linspace(-5,5,100);
x2=linspace(-5,5,100);
[X1,X2]=meshgrid(x1,x2);
f = @(x1,x2) 10*2+(x1.^2-10*cos(2*pi*x1))+(x2.^2-10*cos(2*pi*x2));
f(X1,X2);
Z=f(X1,X2);
contour(X1,X2,Z,0:10:90)
title('Genetic Algorithm','FontWeight','bold','FontSize',20,'FontName','Times New Roman');
xlabel('X_1','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
ylabel('X_2','FontWeight','bold','FontSize',12,'FontName','Times New Roman');
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initial 10 parents
p=20;%number of parents
ns=10;%number of designs
f=@(x) 10*2+(x(1)^2-10*cos(2*pi*x(1)))+(x(2)^2-10*cos(2*pi*x(2)));
X=[0;0];
i=1;
Fminold=80;
c=80;
Sumf=0;%initial sum of fitness
loop=0;
Frank=[0;0];

%initial X
for (i=1:p)
    X(1,i)=rand(1)*10-5;
    X(2,i)=rand(1)*10-5;
        i=i+1;
        
end
plot(X(1,:),X(2,:),'.')%PLOTING
hold on
X%initial X


while loop<200

   Sumf=0; 
      for (i=1:p)
   F(i)=f(X(:,i));
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
    



if loop<=180
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
            ii=p+5;
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
        if f(T(:,i))>f(T(:,j))
            TV=T(:,j);
            T(:,j)=T(:,i);
            T(:,i)=TV;
        end
        j=j+1;
    end
    i=i+1;
end

X=T(:,1:p);

for i=1:p
    
i=i+1;
end



plot(X(1,:),X(2,:),'o')%PLOTING
hold on



loop=loop+1;
end

 Fmin=f(X(:,1))%the min f
 X(:,1)
 plot(X(1,2),X(2,2),'*')
 loop






