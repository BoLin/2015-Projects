%
% the tolerance	for Golden Section			0.001
% following needed for UpperBound_1Var
% the initial guess for UpperBound			lowbound
% the multiplying coeff for UpperBound		a
% the additive coeff for UpperBound	    	b
% The maximum bound for in UpperBound       xmax
% GoldSection_1Var('Example5_1',0.001,0.2,1.6,0.5,10)
%
%  a global statement capturing the values in all iterations is available
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ReturnValue = ...
    GoldSection_1Var(functname,tol,lowbound,a,b,xmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact  % 
warning off  % 
%%% gets upperbound ---------------------------
[bound, fBound] = UpperBound(functname,lowbound,a,b,xmax);


if numel(bound)==3
aL = bound(1);
au = bound(3);
faL = fBound(1);
fau = fBound(3);
else if numel(bound)==2
aL = bound(1);
au = bound(2);
faL = fBound(1);
fau = fBound(2);
    end
end 
%%% sets default tolerance
if (tol == 0)
    tol = 0.0001; %default
end

%%% determines the number of iterations
eps1 = tol/(au - lowbound); %default
tau = 0.38197; % golden ratio
nmax = round(-2.078*log(eps1)) % number of iterations
asL = zeros(nmax,1);
fasL = zeros(nmax,1);
as1 = zeros(nmax,1);
fas1 = zeros(nmax,1);
as2 = zeros(nmax,1);
fas2 = zeros(nmax,1);
asU = zeros(nmax,1);
fasU =zeros(nmax,1);

%%% values at the start iteration
%aL = lowbound;	faL =feval(functname,aL);;
a1 = (1-tau)*aL + tau*au; 
fa1 = feval(functname,a1);

a2 = tau*aL + (1 - tau)*au; 
fa2 = feval(functname,a2);

i = 1;
asL(i) = aL; as1(i) = a1;   as2(i) = a2;  asU(i) = au;
fasL(i) = faL;  fas1(i) = fa1;  fas2(i) = fa2;  fasU(i) = fau;

for i = 1:nmax
    ii = i+1;
    if fa1 >= fa2
        aL = a1;	faL = fa1;
        a1 = a2;	fa1 = fa2;
        a2 = tau*aL + (1 - tau)*au;	fa2 = feval(functname,a2);
        au = au;	fau = fau;  % not necessary -just for clarity
        asL(ii) = aL; as1(ii) = a1;   as2(ii) = a2;  asU(ii) = au;
        fasL(ii) = faL;  fas1(ii) = fa1;  fas2(ii) = fa2;  fasU(ii) = fau;

    else
        au = a2;	fau = fa2;
        a2 = a1;	fa2 = fa1;
        a1 = (1-tau)*aL + tau*au;	fa1 = feval(functname,a1);
        aL = aL;	faL = faL;  % not necessary

        asL(ii) = aL; as1(ii) = a1;   as2(ii) = a2;  asU(ii) = au;
        fasL(ii) = faL;  fas1(ii) = fa1;  fas2(ii) = fa2;  fasU(ii) = fau;

    end
end
% returns the value at the last iteration
%fprintf('xL \t\t FL \t\t x1 \t\t F1 \t\t x2 \t\t F2 \t\t xu \t\t Fu \n');

ReturnValue =[as1(nmax) nmax];