%use subfunction to find objective minimum using 3-point quadratic Polynomial approximation
function [xmin]=quadmin_3point(X1,F1,X2,F2,X3,F3)
%to find the coefficients
a2=(((F3-F1)/(X3-X1))-((F2-F1)/(X2-X1)))/(X3-X2);
a1=((F2-F1)/(X2-X1))-a2*(X1+X2);
a0=F1-a1*X1-a2*X1^2;

%now the minimum for a quadratic approximating polynomial is given by
xmin=-a1/(2*a2);

%fprintf('  \n');
%fprintf('The minimum value is found at %f.\n',xmin);
end