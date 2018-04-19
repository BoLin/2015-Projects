%%Assignment 4 using SQP

function [ ]=SEQ_QUARDHW4(f,g,h,xo,Xl,Xu)
%% using Forward Finite Difference for Sensitivity Calculation
hh=sqrt(2*eps);
gradF=@(x,fvalue) double([(f(x+[hh;0])-fvalue)/hh; (f(x+[0;hh])-fvalue)/hh]);   %using First order forward finite difference
gradg=@(x,gvalue) double([(g(x+[hh;0])-gvalue)/hh; (g(x+[0;hh])-gvalue)/hh]);   %using First order forward finite difference
gradh=@(x,hvalue) double([(h(x+[hh;0])-hvalue)/hh; (h(x+[0;hh])-hvalue)/hh]);   %using First order forward finite difference


%% Now starting the program
x(:,1)=xo ;  %initial values for X!, X2, and X3 slack variable
XXmin=Xl;   XXmax=Xu;
n=2; l=1; jjk=1;
x1(1)=x(1,1);  x2(1)=x(2,1);
funval(1)=f(x(:,1));
consg(1)=g(x(:,1));
consh(1)=h(x(:,1));
B=eye(n) ; Bsub1=B;
it=1;
while(1)
    it=it+1;
    fx=f(x(:,jjk));  delf=double(gradF(x(:,jjk),fx)); deltabar=0.925;
    gx=g(x(:,jjk));  delg=double(gradg(x(:,jjk),gx));     if gx<0; deltaj=1; else deltaj=deltabar; end %here deltaj to check inconsistences in linearized constraints
    hx=h(x(:,jjk));  delh=double(gradh(x(:,jjk),hx));
    Aineq=[delg'];   bineq=[-deltaj*gx];
    Aeq=[delh'];     beq=[-deltabar*hx];
    options = optimset('Algorithm','interior-point-convex');
    S=double(quadprog(B,delf',Aineq,bineq,Aeq,beq,[],[],[],options));

    %direction S found.Now carry out one dimensional seach to minimize phi
    %as an unconstrained function.
    %HERE WE NEED TO FIND THE LANGRANGIAN MULTIPLIERS SO THAT WE CAN USE IT
    %FOR CALCULATION OF AUGMENTED EXTERIOR PENALTY FUNCTION. WE KNOW,
    %L=DELFX+U1*{MAX[0,G(X)]}+U2*ABS(H(X))
    AA=[delg, delf]; bb=[-delf]; lamda(:,jjk)=linsolve(AA,bb);
    if jjk==1
        u1(jjk)=abs(lamda(1)); u2(jjk)=abs(lamda(2));
    elseif jjk>1
        u1prime=u1(jjk-1); u1(jjk)=max(abs(lamda(1)),0.5*(u1prime+abs(lamda(1))));
        u2prime=u2(jjk-1); u2(jjk)=max(abs(lamda(2)),0.5*(u2prime+abs(lamda(2))));
    end
    u11=u1(jjk);   u22=u2(jjk);
    PHI=@(x) f(x)+u11*max(0,g(x))+u22*abs(h(x));
    [alphabst]=optim_alphaHW4(PHI,x(:,jjk),S,XXmin,XXmax);
    if jjk==1
     %   fprintf('\n%d     %f     %f     %f     %f     %f     %f\n',jjk,x(1,jjk),x(2,jjk),f(x(:,jjk)));
    end
    jjk=jjk+1;    %increment in iteration number
    x(:,jjk)=x(:,jjk-1)+alphabst*S ;   %updating the design variables
    x1(it)=x(1,jjk);   x2(it)=x(2,jjk);
    funval(it)=f(x(:,jjk)); consg(it)=g(x(:,jjk)); consh(it)=h(x(:,jjk));
   % fprintf('\n%d     %f     %f     %f     %f     %f     %f\n',jjk,x(1,jjk),x(2,jjk),f(x(:,jjk)));
    
    %% check for convergence.
    %% check for Kuhn-tucker conditions with tolerance
    etol=1e-4;  Gx=g(x(:,jjk));   Hx=h(x(:,jjk));
    if (abs(Gx)<=etol)&&(abs(Hx)<=etol)
        disp('  '); disp('Convergence due to KT conditions satisfied.')
        break
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
    
    if deltaF_absolute>Ea     %to make optimization algorithm robust check the convergence for two or three successive iterations
        N1=0;
    else
        N1=N1+1;    %if satisfied then this N1 is increased
    end
    
    if N1==2;
        disp('  ')
        fprintf('\nConvergence due to absolute change(value=%f) in function criterion in %d iterations\n',deltaF_absolute,jjk)
        break; end
    % return; end     % here the convergence is checked for two successive iterations
    
    %% to check the convergence criteria for absoulute change in objective function
    Er=1e-5;   %tolerance for absolute change in objective function, 0.0001
    deltaF_relative=deltaF_absolute/max(f(x(:,jjk)),10^-10);   %relative change in objective function
    
    if deltaF_relative>Er     %to make optimization algorithm robust check the convergence for two or three successive iterations
        N2=0;
    else
        N2=N2+1;    %if satisfied then this N2 is increased
    end
    
    if N2==2;
        disp('  ')
        fprintf('\nConvergence due to relative change(value=%f) in function criterion in %f iterations\n',deltaF_relative,jjk)
        break; end
    %return; end     % here the convergence is checked for two successive iterations
    
    %% to check the convergence criteria for change in design variables
    vartol=1e-5;
    deltaX=x(:,jjk)-x(:,jjk-1);
    X2norm=norm(deltaX,2);
    if X2norm<=vartol
        % fprintf('\n%d     %f     %f     %f     %f     %f     %f\n',jjk,x(1,jjk),x(2,jjk),f(x(:,jjk)));
        
        disp('  ')
        fprintf('\nConvergence due to change in design variables in %d iterations\n',jjk)
 break; end
    
    
    %% If the convergence criteria are not satisfied then update the matrix B
    fx=f(x(:,jjk));  delf=double(gradF(x(:,jjk),fx));
    gx=g(x(:,jjk)); delg=double(gradg(x(:,jjk),gx));
    hx=h(x(:,jjk));  delh=double(gradh(x(:,jjk),hx));
    AA=[delg, delf];
    bb=[-delf];
    lamda(:,jjk)=linsolve(AA,bb);
    
    p=x(:,jjk)-x(:,jjk-1) ;  %(1)
    
    lamp=@(X,lamda) f(X)+lamda(1)*g(X)+lamda(2)*abs(h(X)) ;   %these values are for new updated x    (2)
    delxlamp=@(X,lamda) [(lamp(X+[hh;0],lamda)-lamp(X,lamda))/hh;(lamp(X+[0;hh],lamda)-lamp(X,lamda))/hh];
    
    
    yy=delxlamp(x(:,jjk),lamda(:,jjk))-delxlamp(x(:,jjk-1),lamda(:,jjk-1))  ;      %(3)
    if p'*yy>=0.2*p'*B*p         %(4)
        theta=1;
    elseif p'*yy<=0.2*p'*B*p
        theta=((0.8*p'*B*p)/(p'*B*p-p'*yy));
    end
    n=theta*yy+(1-theta)*B*p   ;      %(5)
    
    B=B-((B*p*p'*B)/(p'*B*p))+((n*n')/(p'*n)) ; %(6)
    
    if p'*yy>0
        B=B;
    else
        B=Bsub1;
    end
    
end
format short 
mt=[x',funval',consh',consg'];
disp('     x1         x2       f(x)      h(x)      g(x)')
disp(mt);
disp('-----------------------------------------------------------')        
disp('The optimum design obtained from the program are given below:-')
        
        fprintf('\nf(x1,x2)=%f   g(x1,x2)=%f    h(x1,x2)=%f\n',f(x(:,jjk)),g(x(:,jjk)),h(x(:,jjk)));
        fprintf('\nx1=%f   x2=%f\n',x(1,jjk),x(2,jjk));
        
disp('-----------------------------------------------------------')        
        
       

%% %% %% SAVING THE DATA
SQP_f=funval; SQP_g=consg; SQP_h=consh; SQP_xx1=x1;
SQP_xx2=x2; 
save SQP_f.mat; save SQP_g.mat; save SQP_h.mat;  save SQP_xx1.mat; save SQP_xx2.mat;

%% STARTING THE PROGRAM TO PLOT THE FIGURES AND SURFACE PLOT
x1=linspace(0,4,100); x2=linspace(0,4,100);
[X1,X2]=meshgrid(x1,x2);
objecfun=X1.^4-2*X1.^2.*X2+X1.^2+X1.*X2.^2-2.*X1+4;
cons_h= X1.^2+X2.^2-2; cons_g=0.25*X1.^2+.75*X2.^2-1;

%% LOAD THE DATA FOR PLOTTING THE FIGURES
%disp('loading ALM values.......')
load SQP_f.mat ; load SQP_g.mat ;load SQP_h.mat ;load SQP_xx1.mat ; load SQP_xx2.mat ;
%disp('completed!')

%% Plotting the figures WITH CONTOURS AND OPTIMIZATION HISTORY
akkk1=length(SQP_xx1);  akkk2=length(SQP_xx2);
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
contour(X1,X2,cons_h,[SQP_g(8) SQP_g(8)],'b','linewidth',2,'ShowText','on')
contour(X1,X2,cons_g,[SQP_h(8) SQP_h(8)],'g','linewidth',2,'ShowText','on')
plot(SQP_xx1,SQP_xx2,'k','Linewidth',2)
plot(SQP_xx1,SQP_xx2,'-ko','Linewidth',2)

plot(SQP_xx1(akkk1),SQP_xx2(akkk2),'ro','Linewidth',2)
text(SQP_xx1(1)+.03,SQP_xx2(1)+.03, 'xinit','FontName','Times','Fontsize',12,'FontWeight','bold')
text(SQP_xx1(akkk1)+.03,SQP_xx2(akkk2)+.03, 'fmin','FontName','Times','Fontsize',12,'FontWeight','bold')

hold off
end