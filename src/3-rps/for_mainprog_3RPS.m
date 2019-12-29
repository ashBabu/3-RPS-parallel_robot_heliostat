function F=for_mainprog_3RPS(x,OS)
global OP rp zc t Rot

n1=x(1);n3=x(2);o3=x(3);xc=x(4);yc=x(5);

O1G=[xc;yc;zc];   % w.r.t base coordinate system (bcs)
O1P_bcs=Rot'*(OP-t);      % O1P_bcs = O1P described in bcs
g=O1P_bcs-O1G;
GP=g/norm(g,2);
GS=Rot'*OS;         % Sun vector in bcs
A=(GS+GP)/norm((GS+GP),2);
%                 a1=normal(1,1);   a2=normal(2,1);       a3=normal(3,1);      % wrt global coordinate system
n2=-yc/rp;
o1=n2;
o2=n1-(2*xc/rp);

F=[ n1^2    + n2^2    + n3^2    - 1;
    o1^2    + o2^2    + o3^2    - 1;
    A(1)*n1 + A(2)*n2 + A(3)*n3;
    n1*o1   + n2*o2   + n3*o3;
    A(1)*o1 + A(2)*o2 + A(3)*o3
    ];
% F=[ n1^2          + (-yc/rp)^2         + n3^2 - 1;
%     (-yc/rp)^2    + (n1-2*xc/rp)^2     + o3^2 - 1;
%     A(1)*n1       + A(2)*(-yc/rp)      + A(3)*n3;
%     2*n1*(-yc/rp) - 2*(-yc/rp)*(xc/rp) + n3*o3;
%     A(1)*(-yc/rp) + A(2)*(n1-2*xc/rp)  + A(3)*o3;
%    ];


