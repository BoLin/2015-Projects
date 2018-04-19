function [df]=gradientf(functname,X,ii)


h=functname;


    x1=X(1,ii);
    x2=X(2,ii);
    x3=X(3,ii);
    
    f(ii)=h(X(:,ii));
    %x1,x2,x3~=0
            df_dx1 = (h([1.001*x1;x2;x3])-f(ii))/(0.001*x1);
            df_dx2 = (h([x1;1.001*x2;x3])-f(ii))/(0.001*x2);
            df_dx3 = (h([x1;x2;1.001*x3])-f(ii))/(0.001*x3);
    
   df(:,1)=[df_dx1;df_dx2;df_dx3];
   %df(:,1)=[df_dx1;df_dx2];