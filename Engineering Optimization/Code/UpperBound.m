% This function given bounds of minimum of one dimensional problem
% functionname-> objective function in *.m file
% x0-> initial lower guess
% a -> generally 1.61803
% b -> if x0 is zero, this will initiate non zero new value
% bound = UpperBound('Homework2_2',0,1.61803,0.5)
function [bound, fBound] = UpperBound(functname,x0,a,b,xmax)
format compact  % avoid skipping a line when writing to the command window
warning off  % don't report any warnings like divide by zero etc.
%
fl = feval(functname,x0);

if x0==0
    xu = (1+a)*x0+b;
else
    xu = (1+a)*x0;
end
xl = x0;
fu = feval(functname,xu);
if(fu > fl)
       bound = [xl xu];
       fBound = [fl fu];
          disp('Upperbound found but useless');
          disp('modify initial guess');
else
    while (1)
        x1 = xu;
        f1 = fu;
        xu = (1+a)*x1-a*xl;
        if(xu > xmax)
            disp('unbounded');
            break;
        end
        fu=feval(functname,xu);
        if(fu > f1)
            disp('Upperbound found');
            break;
        end
        xl = x1;
        fl = f1;
    end
    bound = [xl x1 xu]
    fBound = [fl f1 fu]
end
