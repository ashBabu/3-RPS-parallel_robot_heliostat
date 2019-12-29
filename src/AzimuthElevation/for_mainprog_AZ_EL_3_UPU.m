function F=for_mainprog_AZ_EL_3_UPU(x,OS)
global OP  OO1 O1M Rot_gam

n1=x(1);n2=x(2);n3=x(3);o1=x(4);o2=x(5); 
o3=0;
O1P_bcs=Rot_gam'*(OP-OO1);      % O1P_bcs = O1P described in bcs
g=O1P_bcs-O1M;
GP=g/norm(g,2);  % unit reflected ray wrt bcs
GS=Rot_gam'*OS;         % Sun vector in bcs
GN=(GS+GP)/norm((GS+GP),2); % unit normal ray in bcs
a1=GN(1); a2=GN(2); a3=GN(3); 

F=[ n1^2    + n2^2    + n3^2    - 1;
    o1^2    + o2^2    + o3^2    - 1;
    a1*n1   + a2*n2   + a3*n3;
    n1*o1   + n2*o2   + n3*o3;
    a1*o1   + a2*o2   + a3*o3
    ];

