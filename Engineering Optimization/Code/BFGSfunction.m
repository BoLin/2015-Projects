function[X]=BFGSfunction(functname,X,ii)

e=1; %initialization:
i=0;
XTEMP(:,1)=X(:,ii);
h=functname;
A=eye(2);
while e>0.001
    
    i=i+1; %counter
    x1=XTEMP(1,i);
    x2=XTEMP(2,i);
    
    f(i)=h(XTEMP(:,i));
    
        if (x1==0) && (x2~=0)
            df_dx1 = (h([0.001;x2])-f(i))/(0.001); 
            df_dx2 = (h([x1;1.001*x2])-f(i))/(0.001*x2);
        else
            df_dx1 = (h([1.001*x1;x2])-f(i))/(0.001*x1);
            df_dx2 = (h([x1;1.001*x2])-f(i))/(0.001*x2);
            
        end
        if (x2==0) && (x1~=0)
            df_dx1 = (h([1.001*x1;x2])-f(i))/(0.001*x1);
            df_dx2 = (h([x1;0.001])-f(i))/(0.001);
        else
            df_dx1 = (h([1.001*x1;x2])-f(i))/(0.001*x1);
            df_dx2 = (h([x1;1.001*x2])-f(i))/(0.001*x2);
        end
        if (x2==0) && (x1==0)
            df_dx1 = (h([0.001;x2])-f(i))/(0.001);
            df_dx2 = (h([x1;0.001])-f(i))/(0.001);
        else
            df_dx1 = (h([1.001*x1;x2])-f(i))/(0.001*x1);
            df_dx2 = (h([x1;1.001*x2])-f(i))/(0.001*x2);
        end
     
    
    df(:,i)=[df_dx1;df_dx2];
     if i==1
        S(:,i)= A*(-df(:,i)/norm(df(:,i)));
    else
        d=XTEMP(:,i)-XTEMP(:,i-1);
        g=df(:,i)-df(:,i-1);
        b=A*d;
        A=A+g*g'/(g'*d)-b*b'/(b'*d);
        S(:,i)= -A^-1*df(:,i);
        S(:,i)=S(:,i)/norm(S(:,i));
     end
     
     [Bracket]=Bracketing_2D_HW3(XTEMP(:,i),S(:,i));
    XTEMP(:,i+1)=Golden_Section_2D_HW3(Bracket,XTEMP(:,i),S(:,i));
    e=max(abs(XTEMP(:,i+1)-XTEMP(:,i)));
     
    if i==60
        break;
    end
    
end
X(:,ii+1)=XTEMP(:,i+1)


