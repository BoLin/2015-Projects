%Assignment 4  
%%AUGMENTED LAGRANGE MULTIPLIER METHOD
function []=AUG_LANGHW4(f,g,h,xo,Xl,Xu)
%f=@(x) x(1,1)^4-2*x(1,1)^2*x(2,1)+x(1,1)^2+x(1,1)*x(2,1)^2-2*x(1,1)+4;
%h=@(x) x(1,1)^2+x(2,1)^2-2;
g1=@(x) g(x);   g2=@(x) -x(1,1);
g3=@(x) x(1,1)-4; g4=@(x) -x(2); g5=@(x) x(2,1)-4;
hh=sqrt(2*eps);
gradF=@(x,fvalue) double([(f(x+[hh;0])-fvalue)./hh; (f(x+[0;hh])-fvalue)./hh]);   %using First order forward finite difference
lamda=ones(1,6);
%% Number of variables and consraints //Initial values
n=2; m=5; l=1;
gamma=3;  %multiplier for rp
rp=1; rp_max=10^5;
X=xo;  %starting value
x1(1)=X(1); x2(1)=X(2);
funval(1)=f(X); consg(1)=g1(X); consh(1)=h(X);
XXmin=Xl;   XXmax=Xu;
Nmax=50;  Ea=1e-4;  Er=1e-4; vartol=1e-4; gradtol=1e-4; tol_constraints=1e-10;
i=1;
disp('iter     x1           x2       f(x1,x2)')
fprintf('\n%d     %6.3f      %6.3f      %6.3f      %6.3f      %6.3f      %6.3f \n',i-1,x1(i),x2(i),funval(i));
while 1
    i=i+1;
    TSI1= max(g1(X),-lamda(1)/(2*rp));
    TSI2= max(g2(X),-lamda(2)/(2*rp));
    TSI3= max(g3(X),-lamda(3)/(2*rp));
    TSI4= max(g4(X),-lamda(4)/(2*rp));
    TSI5= max(g5(X),-lamda(5)/(2*rp));
    TSI=[TSI1,TSI2,TSI3,TSI4,TSI5] ;   %5 TSI's for five inequality constraints
    %%the augmented pseudo-objective function
    A=@(X) f(X)+((TSI(1)*lamda(1)+rp*TSI(1)^2)+(TSI(2)*lamda(2)+rp*TSI(2)^2)+(TSI3*lamda(3)+rp*TSI(3)^2)+(TSI(4)*lamda(4)+rp*TSI(4)^2)+...
        (TSI(5)*lamda(5)+rp*TSI(5)^2))+(lamda(6)*h(X)+rp*h(X)^2);
    X=broyden_goldfarb(A,X,XXmin,XXmax,Nmax,Ea,Er,gradtol);
    %X=broyden_goldfarb(A,X,XXmin,XXmax,Nmax,Ea,Er,gradtol);
    x1(i)=X(1);   x2(i)=X(2);
    funval(i)=f(X);
    consg(i)=g1(X);
    consh(i)=h(X);
    fprintf('\n%d     %6.3f     %6.3f      %6.3f      %6.3f      %6.3f      %6.3f \n',i-1,x1(i),x2(i),funval(i));
    
    %% check for Kuhn-tucker conditions with tolerance
    etol=tol_constraints;  G1=g1(X);  G2=g2(X);   G3=g3(X); G4=g4(X);  G5=g5(X);  Hx=h(X);
    if (abs(G1)<=etol)&&(abs(G2)<=etol)&&(abs(G3)<=etol)&&(abs(G4)<=etol)&&(abs(G5)<=etol)&&(abs(Hx)<=etol)
        
        disp(' '); disp('Convergence due to KT conditions satisfied.')
        break
    end
    
    %% CONVERGENCE CRITERIA
    funcgrad=norm(gradF(X,f(X)));
    if funcgrad<=1e-4
        disp(' ');
        disp('----------------------------------------------------------')
        disp('convergence of grad norm satisfied')
        break ;
    end
    
    %% Convergence due to rp>rpmax
    
    if rp>=rp_max
        disp(' ');
        disp('----------------------------------------------------------')
        disp('convergence satisfied due to rp>rp_max=10^5')
        break
    else
    end
    
    %% updating the langrange multipliers
    %for inequality constraint
    lamda1=lamda(1)+2*rp*max(g1(X),-lamda(1)/(2*rp));
    lamda2=lamda(2)+2*rp*max(g2(X),-lamda(2)/(2*rp));
    lamda3=lamda(3)+2*rp*max(g3(X),-lamda(3)/(2*rp));
    lamda4=lamda(4)+2*rp*max(g4(X),-lamda(4)/(2*rp));
    lamda5=lamda(5)+2*rp*max(g5(X),-lamda(5)/(2*rp));
    lamda6=lamda(6)+2*rp*h(X);
    lamda=[lamda1,lamda2,lamda3,lamda4,lamda5, lamda6];
    
    rp=gamma*rp  ;  %increase rp
    
    
end
disp(' ');
disp('----------------------------------------------------------')
disp('Optimum values:-')
fprintf('\nG1=%f     G2=%f    G3=%f\n',G1,G2,G3);
fprintf('\nG4=%f     G5=%f     Hx=%f\n',G4,G5,Hx);
fprintf('\nx1=%6.3f  x2=%6.3f  f(x)=%6.3f\n',X(1), X(2),f(X));
disp('----------------------------------------------------------')

%% %% %% SAVING THE DATA
ALM_f=funval; ALM_g=consg; ALM_h=consh; ALM_xx1=x1; ALM_xx2=x2;
save ALM_f.mat; save ALM_g.mat; save ALM_h.mat;  save ALM_xx1.mat; save ALM_xx2.mat


%% STARTING THE PROGRAM TO PLOT THE FIGURES AND SURFACE PLOT
x1=linspace(0,4,100);
x2=linspace(0,4,100);
[X1,X2]=meshgrid(x1,x2);
objecfun=X1.^4-2*X1.^2.*X2+X1.^2+X1.*X2.^2-2.*X1+4;
cons_h=  X1.^2+X2.^2-2;
cons_g= 0.25*X1.^2+.75*X2.^2-1;

%% LOAD THE DATA FOR PLOTTING THE FIGURES
disp('loading ALM values.......')
load ALM_f.mat ; load ALM_g.mat ; load ALM_h.mat ; load ALM_xx1.mat ; load ALM_xx2.mat ;
disp('completed!')


%% Plotting the figures WITH CONTOURS AND OPTIMIZATION HISTORY
% figure
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
contour(X1,X2,objecfun,'r','linewidth',1,'ShowText','on')
contour(X1,X2,cons_h,5,'b','linewidth',1,'ShowText','on')
contour(X1,X2,cons_g,5,'g','linewidth',1,'ShowText','on')
legend('f(x)','h(x)','g(x)');
% h2 = ezplot(cons_h,[0,4]);
% set(h2,'Color','b');
% h3 = ezplot(cons_g,[0,4]);
% set(h3,'Color','b');
contour(X1,X2,cons_h,[ALM_h(12),ALM_h(12)],'b','linewidth',2,'ShowText','on')
contour(X1,X2,cons_g,[ALM_g(12),ALM_g(12)],'g','linewidth',2,'ShowText','on')

plot(ALM_xx1,ALM_xx2,'k','Linewidth',2)
plot(ALM_xx1,ALM_xx2,'-ko','Linewidth',2)

akkk1=length(ALM_xx1);  akkk2=length(ALM_xx2);
plot(ALM_xx1(akkk1),ALM_xx2(akkk2),'ro','Linewidth',2)
text(ALM_xx1(1)+.03,ALM_xx2(1)+.03, 'xinit','FontName','Times','Fontsize',12,'FontWeight','bold')
text(ALM_xx1(akkk1)+.03,ALM_xx2(akkk2)+.03, 'fmin','FontName','Times','Fontsize',12,'FontWeight','bold')

hold off
end


%%BFGS METHOD
function [x]=broyden_goldfarb(f,xinitial,XXmin,XXmax,Nmax,Ea,Er,gradtol)

iimax=Nmax;
format long e
%f=@(x) 100*(x(2,1)-x(1,1).^2).^2+100*(1-x(1,1)).^2;

h=sqrt(2*eps);
gradF=@(x,fvalue) double([(f(x+[h;0])-fvalue)./h; (f(x+[0;h])-fvalue)./h]);   %using First order forward finite difference

COUNT=0;   %to count the number of function evaluations

%initial points
x=xinitial;  %zeros(2,1);    %initail values for x1 and x2
fvalue=f(x); COUNT=COUNT+1;
finitial=fvalue;
x1initial=x(1,1);  x2initial=x(2,1);
n=length(x);
objfun_min=zeros(iimax,1);
x1=zeros(iimax,1);
x2=zeros(iimax,1);

%GIving convergence criteria count 
N1=0 ; N2=0;
%fprintf('ii             alpha          x1          x2        f(x1,x2)\n');

%% starting the loop
ii=1;
xqsub1=[x1initial;x2initial]; 
while 1
    if ii==1
        Hq=eye(n);    Hqsub1=Hq;
        delFxqsub1=gradF(x,fvalue);   COUNT=COUNT+2;  %for gradient calculation requires 2 additional functions evaluations for each time
        S=-Hq*delFxqsub1;
        s1(ii)=S(1,1); s2(ii)=S(2,1);
        [alphaoptim]=optim_alphaHW4(f,x,S,XXmin,XXmax);    %THis subalgorithm finds the optimum value for alpha during each iteration (alphamin=0,alphamax=5)
                x=x+alphaoptim*[s1(ii);s2(ii)];
       % COUNT=COUNT+fun_countalpha;   %COUNTING THE FUNCTION EVALUATIONS FOR ONE DIMENSIONAL SEARCH
        if alphaoptim==0; break; end %return; end    
        x1(ii)=x(1,1);
        x2(ii)=x(2,1);

         x=[x1(ii);x2(ii)];
        fvalue=f(x);            COUNT=COUNT+1;     %here one more CALCULATION
        objfun_min(ii)=fvalue;
        
        delFxq=gradF(x,fvalue);      COUNT=COUNT+2;  %for gradient calculation requires 2 additional functions evaluations
        xq=[x1(ii);x2(ii)]; %for ii==1 xq-1 is 0.. computer doesnt calculate index zero so use this
        DELTAX=xq-xqsub1;
       YK=delFxq-delFxqsub1;
        Dq=symD(xq,xqsub1,delFxq,delFxqsub1,Hq);
        Hq=Hq+Dq;
        
        %CHECKING POSITIVE DEFINITENESS
        MTP=eig(Hqsub1); 
      if (MTP(1,1)>0)&&(MTP(2,1)>0)&&(DELTAX'*YK>0)
          Hq=Hq;
      else
          Hq=Hqsub1;
      end
        
        
        delFxqsub1=delFxq;
       % fprintf('\n%f     %f     %f     %f     %f     %f     %f\n',ii,alphaoptim,x1(ii),x2(ii),objfun_min(ii));
        
    else
        S=-Hq*delFxqsub1;     Hqsub1=Hq;
        s1(ii)=S(1,1); s2(ii)=S(2,1);     
           [alphaoptim]=optim_alphaHW4(f,x,S,XXmin,XXmax);    %THis subalgorithm finds the optimum value for alpha during each iteration (alphamin=0,alphamax=5)
                x=x+alphaoptim*[s1(ii);s2(ii)];
       % COUNT=COUNT+fun_countalpha;   %COUNTING THE FUNCTION EVALUATIONS FOR ONE DIMENSIONAL SEARCH                               
        if alphaoptim==0; break; end %return; end
        x1(ii)=x(1,1);
        x2(ii)=x(2,1);
       
         x=[x1(ii);x2(ii)];
        fvalue=f(x);            COUNT=COUNT+1;     %here one more CALCULATION
        objfun_min(ii)=fvalue;
        delFxq=gradF(x,fvalue);   COUNT=COUNT+2;  %for gradient calculation requires 2 additional functions evaluations
        xq=[x1(ii);x2(ii)]; xqsub1=[x1(ii-1);x2(ii-1)];  %for ii==1 xq-1 is 0.. computer doesnt calculate index zero so use this
        DELTAX=xq-xqsub1;
       YK=delFxq-delFxqsub1;
        Dq=symD(xq,xqsub1,delFxq,delFxqsub1,Hq);
        Hq=Hq+Dq;
        
                %CHECKING POSITIVE DEFINITENESS
        MTP=eig(Hqsub1); 
      if (MTP(1,1)>0)&&(MTP(2,1)>0)&&(DELTAX'*YK>0)
          Hq=Hq;
      else
          Hq=Hqsub1;
      end
        
        
        delFxqsub1=delFxq;
      %  fprintf('\n%f     %f     %f     %f     %f     %f     %f\n',ii,alphaoptim,x1(ii),x2(ii),objfun_min(ii));
        
        
        
         %% to check the convergence criteria for absoulute change in objective function
        % Ea=1e-10;   %tolerance for absolute change in objective function, 0.0001
        
        FXq=objfun_min(ii); FXq_1=objfun_min(ii-1);
        
        deltaF_absolute=abs(FXq-FXq_1);
        
        if deltaF_absolute>Ea     %to make optimization algorithm robust check the convergence for two or three successive iterations
            N1=0;
        else
            N1=N1+1;    %if satisfied then this N1 is increased
        end
        if N1==2;
        %    fprintf('\nConvergence due to absolute change in function criterion in %d iterations\n',ii)
            break; end
        %return; end     % here the convergence is checked for two successive iterations
        
        
        %% to check the convergence criteria for absoulute change in objective function
        % Er=1e-10;   %tolerance for relative change in objective function, 0.0001
        deltaF_relative=deltaF_absolute/max(FXq,10^-10);   %relative change in objective function
        
        if deltaF_relative>Er     %to make optimization algorithm robust check the convergence for two or three successive iterations
            N2=0;
        else
            N2=N2+1;    %if satisfied then this N2 is increased
        end
        if N2==2;
       %     fprintf('\nConvergence due to relative change in function criterion in %d iterations\n',ii)
            break; end
        %return; end     % here the convergence is checked for two successive iterations
        
        
    end
    
    
            %%     %%    %% to check the convergence criteria for norm of  gradient
        
        gradFnorm=norm(delFxq,2);
        if gradFnorm<=gradtol
           % fprintf('convergence criteria satisfied due to norm of gradient of function(value=%f)\n',gradFnorm)
            break; end
        
    if ii>=iimax
      %  disp('maximum number of iterations exceeded')
        break;
    end
    ii=ii+1;
end


end

%find the minimum of univariate function using GOLDEN-SECTION METHOD
function [y]=optim_alphaHW4(f,xo,S,Xmin,Xmax); 
%function [y,fun_countalpha]=optim_alpha(f,xo,S,Xmin,Xmax); 
xo;
S;
if sign(S(1))==1,
  x1u=(Xmax(1)-xo(1))/S(1); %S(1)>0
else
    x1u=(Xmin(1)-xo(1))/S(1)  ;    %S(1)<0 
end

if sign(S(2))==1, 
    x2u=(Xmax(2)-xo(2))/S(2);   %S(2)>0
else
    x2u=(Xmin(1)-xo(2))/S(2); 
end          %S(2)><0
    
ALPHA=[x1u,x2u];

alphamax=min(ALPHA(find(ALPHA>0)));

alphamin=0;

fun_countgs=0;
fun_countalpha=0;

[Al,Fl,Au,Fu,fun_countbrkt]=bounds(f,xo,S,alphamin,alphamax);
%we can also three points Xl,X1,Xu found from the bracketing of a minima.
eps=1e-3;   %supposing the final interval is 0.1% of the initial
tau=0.381966;
N=round(log(eps)/log(1-tau));
Aone=(1-tau)*Al+tau*Au;             
Fone=f(xo+Aone*S);       
Atwo=tau*Al+(1-tau)*Au;
Ftwo=f(xo+Atwo*S);

fun_countalpha=fun_countbrkt+2;
K=0;
fun_countgs=2;
while 1
    K=K+1;
    if Fone>Ftwo
        Al=Aone; Fl=Fone;
        Aone=Atwo; Fone=Ftwo;
        Atwo=tau*Al+(1-tau)*Au; Ftwo=f(xo+Atwo*S); fun_countgs=fun_countgs+1; fun_countalpha=fun_countalpha+1;
    elseif Fone<Ftwo
        Au=Atwo; Fu=Ftwo;
        Atwo=Aone; Ftwo=Fone;
        Aone=(1-tau)*Al+tau*Au;   Fone=f(xo+Aone*S); fun_countalpha=fun_countalpha+1; fun_countgs=fun_countgs+1;
        
    end
    %  fprintf('%f     %f     %f     %f     %f     %f\n',Al,Fl,Aone,Fone,Atwo,Ftwo,Au,Fu);
    if K>=N
        %   fprintf('maximum number of iterations %d exceeded.-->StOPPED',N);
        if Fone>Ftwo
            Amin=Atwo; Fmin=Ftwo;
        else
            Amin=Aone; Fmin=Fone;
        end
        %fprintf('\n      \n');
        % fprintf('The best approximation to xmin and fmin using Golden-Section method are:- \n');
        % fprintf('Amin=%f       Fmin=%f\n',Amin,Fmin);
        break;
    end
end
y=Amin;
Fmin;
K;

fun_countgs;
fun_countalpha;
%we can again refine this solution using 3 point quadratic or 4 point cubic
%xmin=cubic_4point(Xl,Fl,Xone,Fone,Xtwo,Ftwo,Xu,Fu);
end



function [Al,Fl,Au,Fu,fun_countbrkt]=bounds(f,xo,S,alphamin,alphamax);   %here chosing a proper interval for alpha

a=1.61803;  %(Golden section ratio)

Al=alphamin;
fun_countbrkt=0;
Fl=f(xo+Al*S);   %  dXl1= dX1(Xl);   dXl2=dX2(Xl);
S;
Amax=alphamax;
deltaA=.1*Amax; %abs(.3*f(Xl)/dXl1);      %here increases the step size higher
Au=Al+deltaA; %Xl+deltax; if we know its unimodal and decreasing
Fu=f(xo+Au*S);
fun_countbrkt=2;
while 1    %use the upper limit Xmax
    if Fu>Fl
        break;
    end
    Aone=Au;
    Fone=Fu;
    Au=(1+a)*Aone-a*Al;
    if Au>Amax
       % Au=Amax
       %Fu=f(xo+Au*S)
       %disp('unbounded')
        break;
    end
    Fu=f(xo+Au*S);
    fun_countbrkt=fun_countbrkt+1;
    if Fu>Fone
        break;
    end
    Al=Aone;
    Fl=Fone;
    
    %if Al>=alphamax
      %  break
   % end
end
fun_countbrkt;

end





%%   %to caluclate D for BFGS
function Dq=symD(xq,xqsub1,delFxq,delFxqsub1,Hq)
p=xq-xqsub1;
y=delFxq-delFxqsub1;
sigma=p'*y;
tau=y'*Hq*y;
Dq=((sigma+tau)/(sigma^2))*p*p'-(1/sigma)*(Hq*y*p'+p*(Hq*y)');
end

