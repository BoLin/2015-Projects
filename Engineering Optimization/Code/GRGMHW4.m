%%Assignment 4
%%General Reduced Gradient Method
function [ ]=GRGM(f,g,h,xo,Xl,Xu)
%now the variable x=[X1;X2;X3]
h1=@(x) x(1)^2+x(2)^2-2;   %equality constraint for main problem
h2=@(x) g([x(1);x(2)])+x(3);  %creating a equality constraint from inequality constraint g
%the side constraints are given as :-
%0<=X1<=4
%0<=X2<=4
%0<=X3<=1e5   % upper limit for slack variables is set very large(infinity)
%% using Forward Finite Difference for Sensitivity Calculation
hh=sqrt(2*eps);
gradF=@(x,fvalue) double([(f(x+[hh;0;0])-fvalue)./hh; (f(x+[0;hh;0])-fvalue)./hh; (f(x+[0;0;hh])-fvalue)./hh]);   %using First order forward finite difference
gradh1=@(x,h1value) double([(h1(x+[hh;0;0])-h1value)./hh; (h1(x+[0;hh;0])-h1value)./hh; (h1(x+[0;0;hh])-h1value)./hh]);   %using First order forward finite difference
gradh2=@(x,h2value) double([(h2(x+[hh;0;0])-h2value)./hh; (h2(x+[0;hh;0])-h2value)./hh;  (h2(x+[0;0;hh])-h2value)./hh]);   %using First order forward finite difference

%% Now starting the program
n=2; l=1; m=1;
x=xo;
x1(1)=x(1,1);  x2(1)=x(2,1);
funval(1)=f(x(:,1));
consg(1)=g(x(:,1));
consh(1)=h(x(:,1));
h1value=h1(x(:,1));
h2value=h2(x(:,1));
depvar=[2 3]; Y=x(depvar,1)  ;  %dependent variables
varnum=[1:n+m]; varnum(depvar)=[]; indvar=varnum; 
Z=x(indvar,1) ;%independant variables
fprintf('\n The independent variables are X%d.\n',indvar)
fprintf('\n The dependent variables are X%d.\n',depvar)
it=1; jjk=1;
%% STARTING THE LOOP
while 1
    it=it+1;
    %Calculate the gradients of objective functions and all constraints
    fvalue=f(x(:,jjk)); h1value=h1(x(:,jjk)); h2value=h2(x(:,jjk));
    DELFX=gradF(x(:,jjk),fvalue); DELh1X=gradh1(x(:,jjk),h1value); DELh2X=gradh2(x(:,jjk),h2value);
    
    %% finding the gradients of independant and dependant variables separately
    delZ_FX=DELFX(indvar);   delY_FX=DELFX(depvar);
    delZ_h1X=DELh1X(indvar); delY_h1X=DELh1X(depvar);
    delZ_h2X=DELh2X(indvar); delY_h2X=DELh2X(depvar);
    A=[delZ_h1X';delZ_h2X'];
    B=[delY_h1X';delY_h2X'];
    GR=delZ_FX-(inv(B)*A)'*delY_FX ;   %checking the GR from QX=I and B^-1*A
    %% search direction either can use S=-GR or can use CGM
    S=-GR  ;   %for independant variable
    for jjjk=1:length(S)
        if sign(S(jjjk))==1
            alphaZ(jjjk)=(Xu(indvar(jjjk))-Z(jjjk))/S(jjjk)  ;  %gives toward upper bound for independant variables
        else
            alphaZ(jjjk)=(Xl(indvar(jjjk))-Z(jjjk))/S(jjjk) ;   %gives toward lower bound  for independant variables
        end
        jjjk=jjjk+1;
    end
    %for dependant variables
    Yp=-inv(B)*A*S;
    for jjjk=1:length(Yp)
        if sign(Yp(jjjk))==1
            alphaY(jjjk)=(Xu(depvar(jjjk))-Y(jjjk))/Yp(jjjk) ;  %gives toward upper bound for dependant variables
        else
            alphaY(jjjk)=(Xl(depvar(jjjk))-Y(jjjk))/Yp(jjjk) ;   %gives toward lower bound  for dependant variables
        end
        jjjk=jjjk+1;
    end
    
    alpha(1)=0;
    %alpha(3)=min([alphaZ alphaY])   %gives the suitable value for step parameter
    alpha(3)=.063; alpha(2)=(alpha(1)+alpha(3))/2; 
    XXxx=x(:,jjk);
    
    %% DELETE SOME VALUES HERE USED FOR CHECKING ONLY
    for tt=1:3
        X=XXxx;  X(indvar)=X(indvar)+alpha(tt)*S;
        dY=-inv(B)*(A*alpha(tt)*S) ;
        mmt=1; xtst(:,mmt)=X;
        while(1)
            mmt=mmt+1;   
            X(depvar)=X(depvar)+dY;
            if X(1)<=0; X(1)=0;  end
            if X(2)<=0; X(2)=0; end
            if X(3)<=0; % X(3)=0;
            end
            xtst(:,mmt)=X; FFX=f(X); h1x=h1(X); h2x=h2(X);
            if (h1x<=1e-4)&&(h2x<=1e-4)
                % disp('OK NEWTON METHOD!!')
                alpha(tt)=alpha(tt); FFXX(tt)=FFX; tt=tt+1;
                break
            else
                dY=inv(B)*(-[h1(X);h2(X)])    ; 
            end
        end
    end
    [alphabst]=quadmin_3pointHW4(alpha(1),FFXX(1),alpha(2),FFXX(2),alpha(3),FFXX(3));
    jjk=jjk+1;        %increasing the  iteration number
    x(indvar,jjk)=x(indvar,jjk-1)+alphabst*S;
    x(depvar,jjk)=X(depvar);
    %dY=inv(B)*(-[h1(x(:,jjk));h2(x(:,jjk))]-A*alphabst*S) ;
    dY=-inv(B)*(A*alphabst*S) ;
    x(depvar,jjk)=x(depvar,jjk)+dY;
    if x(1,jjk)<0; x(1,jjk)=0; end
    if x(2,jjk)<0; x(2,jjk)=0; end
    x1(it)=x(1,jjk);   x2(it)=x(2,jjk);
    funval(it)=f(x(:,jjk)); consg(it)=g(x(:,jjk)); consh(it)=h(x(:,jjk));
    %% check for convergence.
    %% check for Kuhn-tucker conditions with tolerance
    etol=1e-3;  Gx=g(x(:,jjk));   Hx=h(x(:,jjk));
    if (abs(Gx)<=etol)&&(abs(Hx)<=etol)
        disp('  '); disp('Convergence due to KT conditions satisfied.')
        break;
    end
    
    %% CONVERGENCE CRITERIA
    funcgrad=norm(gradF(x(:,jjk),f(x(:,jjk))));
    if funcgrad<=1e-4
        disp('  '); disp('convergence of grad norm satisfied')
        break ;
    end
    
    %% to check the convergence criteria for absoulute change in objective function
    Ea=1e-4;   %tolerance for absolute change in objective function, 0.0001
    deltaF_absolute=abs(f(x(:,jjk))-f(x(:,jjk-1)));
    if deltaF_absolute>Ea; N1=0; else N1=N1+1; end
    if N1==2; disp('  ')
        fprintf('\nConvergence due to absolute change(value=%d) in function criterion in %d iterations\n',deltaF_absolute,jjk);
        break; end
    
    %% to check the convergence criteria for absoulute change in objective function
    Er=1e-4;   %tolerance for absolute change in objective function, 0.0001
    deltaF_relative=deltaF_absolute/max(f(x(:,jjk)),10^-10);   %relative change in objective function
    if deltaF_relative>Er; N2=0; else N2=N2+1; end
    if N2==2;
        disp('  ')
        fprintf('\nConvergence due to relative change(value=%f) in function criterion in %d iterations\n',deltaF_relative,jjk)
        break;
    end
    
    %% to check the convergence criteria for change in design variables
    vartol=1e-4;
    deltaX=x(:,jjk)-x(:,jjk-1); X2norm=norm(deltaX,2);
    if X2norm<=vartol
        disp('  ')
        fprintf('\nConvergence due to change in design variables in %d iterations\n',jjk)
        break;
    end
end

disp('----------------------------------------------------------------')
disp('The optimum design obtained from the program are given below:-')

fprintf('\nObjective function: f(x1,x2)=%f\n',f(x(:,jjk)));

fprintf('\nConstraints:        g(x1,x2)=%f    h(x1,x2)=%f\n',g(x(:,jjk)),h(x(:,jjk)));


fprintf('\nOptimum design values: x1=%f   x2=%f\n',x(1,jjk),x(2,jjk));
disp('----------------------------------------------------------------')

format short


mt=[x',funval',consh',consg'];
disp('     x1         x2        x3       f(x)      h(x)      g(x)')
disp(mt);

%% %% %% SAVING THE DATA
GRGM_f=funval; GRGM_g=consg; GRGM_h=consh; GRGM_xx1=x1; GRGM_xx2=x2;
save GRGM_f.mat ; save GRGM_g.mat; save GRGM_h.mat; save GRGM_xx1.mat; save GRGM_xx2.mat

%% STARTING THE PROGRAM TO PLOT THE FIGURES AND SURFACE PLOT
x1=linspace(0,1.5,100); x2=linspace(0,3,100);
%x1=linspace(0,2,100); x2=linspace(0,3,100);   %chose this for another
%design variable
[X1,X2]=meshgrid(x1,x2);
objecfun=X1.^4-2*X1.^2.*X2+X1.^2+X1.*X2.^2-2.*X1+4;
cons_h= X1.^2+X2.^2-2; cons_g=0.25*X1.^2+.75*X2.^2-1;
%% LOAD THE DATA FOR PLOTTING THE FIGURES
disp('loading GRGM values.......')
load GRGM_f.mat ; load GRGM_g.mat ;load GRGM_h.mat ;load GRGM_xx1.mat ; load GRGM_xx2.mat ;
disp('completed!')

%% Plotting the figures WITH CONTOURS AND OPTIMIZATION HISTORY
akkk1=length(GRGM_xx1);  akkk2=length(GRGM_xx2);
figure
ax=gca;
set(ax,'FontName','Times','Fontsize',12,'FontWeight','bold','Fontangle','italic');
box on
grid
set(gcf,'color',[1 1 1])   %to make the backgroung white
title('Optimization history of function f(x)=x1^4-2*x1^2*x2+x1^2+x1*x2^2-2*x1+4',...
    'FontName','Times','Fontsize',14,'FontWeight','bold');
xlabel('x1','FontName','Times','Fontsize',14,'FontWeight','bold');
ylabel('x2','FontName','Times','Fontsize',14,'FontWeight','bold');

hold on
contour(X1,X2,objecfun,5,'r','linewidth',1,'ShowText','on')
contour(X1,X2,cons_h,5,'b','linewidth',1,'ShowText','on')
contour(X1,X2,cons_g,5,'g','linewidth',1,'ShowText','on')

legend('f(x)','h(x)','g(x)');
contour(X1,X2,objecfun,[GRGM_f(akkk1) GRGM_f(akkk1)],'r','linewidth',2)
contour(X1,X2,cons_h,[GRGM_h(akkk1) GRGM_h(akkk1)],'b','linewidth',2)
contour(X1,X2,cons_g,[GRGM_g(akkk1) GRGM_g(akkk1)],'g','linewidth',2)
plot(GRGM_xx1,GRGM_xx2,'k','Linewidth',2)
plot(GRGM_xx1,GRGM_xx2,'-ko','Linewidth',2)
plot(GRGM_xx1(akkk1),GRGM_xx2(akkk2),'ro','Linewidth',2)
text(GRGM_xx1(1)+.03,GRGM_xx2(1)+.03, 'xinit','FontName','Times','Fontsize',12,'FontWeight','bold')
text(GRGM_xx1(akkk1)+.03,GRGM_xx2(akkk2)+.03, 'fmin','FontName','Times','Fontsize',12,'FontWeight','bold')
hold off
end