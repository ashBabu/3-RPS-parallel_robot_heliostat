function [c,ceq] = confuneq_error_min(X)
global rp 
n1 = X(1); n3=X(2); o3 = X(3); a1N=X(4);   a2N = X(5);   a3N = X(6);
xc = X(7); yc = X(8);
n2 = -yc/rp;
o1 = n2;
o2 = n1 - (2*xc/rp);
% Nonlinear inequality constraints
c = [];
% Nonlinear equality constraints
% ceq = [n1^2  + n2^2  + n3^2 -1;
%        o1^2  + o2^2  + o3^2 -1;
%        a1N^2 + a2N^2 + a3N^2-1;
%        n1*o1 + n2*o2 + n3*o3;
%        a1N*n1 + a2N*n2 + a3N*n3;
%        a1N*o1 + a2N*o2 + a3N*o3];

       ceq(1) = n1^2  + n2^2  + n3^2 -1;
       ceq(2) = o1^2  + o2^2  + o3^2 -1;
       ceq(3) = a1N^2 + a2N^2 + a3N^2-1;
       ceq(4) = n1*o1 + n2*o2 + n3*o3;
       ceq(5) = a1N*n1 + a2N*n2 + a3N*n3;
       ceq(6) = a1N*o1 + a2N*o2 + a3N*o3;
%        ceq(7) = n2 + yc/rp;
%        ceq(8) = o1- n2;
%        ceq(9) = o2 - n1+ (2*xc/rp);
       