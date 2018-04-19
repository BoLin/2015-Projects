function f=pred(x)
global ModelInfo
X=ModelInfo.X;
y=ModelInfo.y;
theta=10.^ModelInfo.Theta;
U=ModelInfo.U;
n=size(X,1);
one=ones(n,1);
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
psi=ones(n,1);
for i=1:n
    psi(i)=exp(-sum(theta.*abs(X(i,:)-x).^2));
end
f=mu+psi'*(U\(U'\(y-one*mu)));