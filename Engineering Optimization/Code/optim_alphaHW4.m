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