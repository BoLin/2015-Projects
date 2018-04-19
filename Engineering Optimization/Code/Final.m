
global ModelInfo
ModelInfo.X=[0;3;4;6];
f=@(x) sin(0.1+2*x)/(0.1+x)
Y=[f(0);f(3);f(4);f(6)]
UpperTheta=ones(1,1).*2;
LowerTheta=ones(1,1).*-3;
ModelInfo.y(:,1)=Y;
[ModelInfo.Theta,Min]=ga(@likelihood,1,[],[],[],[],LowerTheta,UpperTheta);
[NegLnLike,ModelInfo.Psi,ModelInfo.U]=likelihood(ModelInfo.Theta);

[x,fmin]=ga(@pred,1,[],[],[],[],0,6)
S=std2([Y;fmin])
A=(f(6)-fmin)/S
a=cdf('norm',A,f(6))
b=pdf('norm',A,mean([Y;fmin]),S)
E=(f(6)-fmin)*1+S*b

f=@(x) sin(0.1+2*x)/(0.1+x)
X=[0;2;3;4;5;6]
Y=[f(0);f(2);f(3);f(4);f(5);f(6)]


X=[1 0 0 0;1 2 4 8;1 3 9 27;1 4 16 64;1 5 25 125;1 6 36 216]
B=inv(X'*X)*(X'*Y)
SS=Y'*Y-B'*X'*Y
rms=sqrt(SS/5)
fhat=@(x) B'*[1;x;x.^2;x.^3]
Yhat=[fhat(0);fhat(3);fhat(4);fhat(5);fhat(6)]
Standarderror=std(Yhat)
e=(Y-Yhat)
Sumofsquareerror=e'*e

Sy=(Y-(ones(1,5)*Y/5)*ones(5,1))
SSy=Sy'*Sy

Sr=(Yhat-(ones(1,5)*Y/5)*ones(5,1))
SSr=Sr'*Sr
Rsq=SSr/SSy
Ra=1-(1-Rsq^2)*(4/1)

ny=5;nbeta=4;
bnds=[0;6]
[dce,x]=cordexch(1,6,'q','bounds',bnds);

X=[1 0 0 0;1 3 9 27;1 4 16 64;1 5 25 125;1 6 36 216;1 1 1 1 ]
