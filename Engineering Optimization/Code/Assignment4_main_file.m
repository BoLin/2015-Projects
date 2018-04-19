clc
clear all
f=@(x) x(1,1)^4-2*x(1,1)^2*x(2,1)+x(1,1)^2+x(1,1)*x(2,1)^2-2*x(1,1)+4;
g=@(x) 0.25*x(1,1)^2+.75*x(2,1)^2-1;   %inequality constraint
h=@(x) x(1,1)^2+x(2,1)^2-2;   %original equlaity constraint
Xl=[0;0];
Xu=[4;4];
xo=[3;2];     % CHOOSING INITIAL DESIGN VECTOR FOR OPTIMIZATION
%% IN ORDER TO RUN ALM, SQP, AND GRGM. PLEASE REMOVE THE "%" AND RUN THE PROGRAM INDIVIDUALLY


%AUGMENTED LANGRANGIAN METHOD
%AUG_LANGHW4(f,g,h,xo,Xl,Xu)

%SEQUENTIAL QUADRATIC PROGRAMMING 
SEQ_QUARDHW4(f,g,h,xo,Xl,Xu)

%GENERALLY REDUCED GRADENT METHOD (THREE DESIGN POINTS GIVEN)
% xo=[1.4;2;.3450];
% xo=[3;2;.3450];         %this does not provide optimal solution so another
% %design point is chosen as provided in the  Book
%  xo=[1.3;.5568;.3450];
%  Xl=[0;0;0]; Xu=[4;4;1e5];   % USE THIS LIMIT IN ADDTION FOR SLACK VARIABLES
%  GRGMHW4(f,g,h,xo,Xl,Xu)