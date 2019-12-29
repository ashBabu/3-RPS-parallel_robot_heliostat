function distN = error_min(X,a,b,c,OS,zeta)
    global rp t Rot m  zc
    n1N = X(1); n3N = X(2); o3N = X(3);
    a1N = X(4); a2N = X(5); a3N = X(6);      % wrt base coordinate system
    xcN = X(7); ycN = X(8);
    n2N = -ycN/rp;
    o1N = n2N;
    o2N = n1N - (2*xcN/rp);
    O1GN = [xcN;ycN;zc];  % new G in bcs
    OGN = t+Rot*O1GN;   % This is the new point G in gcs
    RN  = [n1N o1N a1N ;               % Error corrected rotation matrix which 
           n2N o2N a2N ;               % takes mirror to base coordinate system
           n3N o3N a3N ];
    
    Gm1= [m;m;0]; Gm2= [-m;m;0]; Gm3= [-m;-m;0]; Gm4= [m;-m;0];  % G is the centre of the mirror.
    Om1N = OGN + RN*Gm1; Om2N = OGN + RN*Gm2;  %% mirror corners in gcs
    Om3N = OGN + RN*Gm3; Om4N = OGN + RN*Gm4;
   
    %%  New reflected ray
    axy = cross(Rot'*OS,[a1N;a2N;a3N]);  % IN bcs
    axs = axy/norm(axy,2);   % IN bcs
    angle1 = acos(dot(Rot'*OS,[a1N;a2N;a3N]));  % IN bcs);  %%% RADIANS
    k1=axs(1); k2 = axs(2); k3 = axs(3);
    
    r11 = k1^2*(1-cos(angle1)) +  cos(angle1);
    r12 = k1*k2*(1-cos(angle1))-  k3*sin(angle1);              
    r13 = k3*k1*(1-cos(angle1))+  k2*sin(angle1);
    r21 = k1*k2*(1-cos(angle1))+  k3*sin(angle1);
    r22 = k2^2*(1-cos(angle1)) +  cos(angle1);                  
    r23 = k3*k2*(1-cos(angle1))-  k1*sin(angle1);
    r31 = k3*k1*(1-cos(angle1))-  k2*sin(angle1);
    r32 = k3*k2*(1-cos(angle1))+  k1*sin(angle1);
    r33 = k3^2*(1-cos(angle1)) +  cos(angle1);
    Rot_ref = [r11 r12 r13;
               r21 r22 r23;
               r31 r32 r33];
    GPN = Rot_ref*[a1N;a2N;a3N] ;
    
%      fhandle=@(r)Refl_ray(r,Rot'*OS,[a1N;a2N;a3N]);
%     [GPN,fvald] = fsolve(fhandle,r0);
% %     r0=GPN;
%     %%%%
    %%
    
    k1N=(cosd(zeta)*(a-Om1N(1)) + sind(zeta)* (b-Om1N(2)))/( GPN(1)*cosd(zeta) + GPN(2)*sind(zeta) );
    k2N=(cosd(zeta)*(a-Om2N(1)) + sind(zeta)* (b-Om2N(2)))/( GPN(1)*cosd(zeta) + GPN(2)*sind(zeta) );
    k3N=(cosd(zeta)*(a-Om3N(1)) + sind(zeta)* (b-Om3N(2)))/( GPN(1)*cosd(zeta) + GPN(2)*sind(zeta) );
    k4N=(cosd(zeta)*(a-Om4N(1)) + sind(zeta)* (b-Om4N(2)))/( GPN(1)*cosd(zeta) + GPN(2)*sind(zeta) );
   
    % hitting points on the receiver
    hp1N=Om1N+k1N*GPN; % in gcs
    hp2N=Om2N+k2N*GPN; % in gcs
    hp3N=Om3N+k3N*GPN; % in gcs
    hp4N=Om4N+k4N*GPN; % in gcs
    centroidN = (hp1N+hp2N+hp3N+hp4N)/4;
    distN = norm((centroidN-[a;b;c]),2);
end