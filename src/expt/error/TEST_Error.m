close all 
clear all
clc
home
set(0,'defaulttextinterpreter','latex')
% set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 25)
set(0,'DefaultAxesFontWeight', 'Bold')

% Change default text fonts.
% set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 13)

%%
% % Bangalore Latitude
lad=12;
lam=58;
las=13;

% Bangalore Longitude
lod=77;
lom=33;
los=37;

latitude=lad+(lam+las/60)/60;  % in degrees
longitude=lod+(lom+los/60)/60; % in degrees
%%
% Rajasthan Coordinates
% latitude = 27.0238;
% longitude = 74.2179;
%%

year=2016;
month=3;
day=20;
d=date(year,month,day);

time_incre=360;
LT=0:time_incre:24*3600;

global OP rp zc  Rot t m
psi= 0;
rad= 50;
tx=rad*cosd(psi);ty=rad*sind(psi);tz=0; t=[tx;ty;tz];
a=0;b=0;c=25;
% OP=[a;b;c];
%%
zeta=psi; % the normal (lying in a plane parallel to XY plane) to the receiver is at an angle zeta to the X axis 
% Trans_Rec = Rz_90*Rx_90*Ry_zeta

Rz_90= [ cosd(90)  -sind(90) 0 ;
         sind(90)   cosd(90) 0;
            0          0     1];
Rx_90 =   [1       0          0    ;
           0   cosd(90)   -sind(90);
           0   sind(90)   cosd(90)];
Ry_zeta =   [cosd(zeta) 0  sind(zeta);
                0       1      0     ;
            -sind(zeta) 0  cosd(zeta)];      
       
Rot_Rec = Rz_90*Rx_90*Ry_zeta;
Trans_Rec = [Rot_Rec,[a;b;c];[0 0 0 1]];
%%

OP_rec_coord = [0.15;0.35;0.0];  % OP in receiver coordinates.
OP = Trans_Rec * [OP_rec_coord;1]; % OP in gcs
OP(4)=[];
zc=2;
% zc = 0.2365;
alpha=120;
beta=120;
mirror_size = 1.6; % 2m x 2m square mirror
rp = mirror_size/4;  
rb = rp;
[elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude);
%% Plot of azimuth and elevation angles
v=8/time_incre*3600+1; % start time for the plots
w=17/time_incre*3600+1; % end time for the plots

figure(1)
plot(LT(v:w)/3600,elevation((v: w)),'r','LineWidth',4)
hold on
plot(LT(v:w)/3600,azimuth((v: w)),'b','LineWidth',4)
xlabel('Time  \it{t / hr}')
ylabel('$Sun vector angles ~ \it{\theta / (^0)}$')
lgnd=legend('Elevation',' Azimuth');
grid on
%%
[xs ,ys, zs]=sunposition(elevation, azimuth);  % xs,ys,zs are coordinates of the sun
OS=[xs;ys;zs];

gamma=0;       % Rotation of the base platform w.r.t the global coordinate system.
Rot=[ cosd(gamma) -sind(gamma) 0;
      sind(gamma)  cosd(gamma) 0;
          0            0       1]; 
x0 = [.7;-.6;0.1;-0.7;-0.01]; % works for psi=0,45,90,135,180,225,270,315 degrees

% x0 = [0.2;0.61;0.5;-0.01;0.015];  % works only for psi=45 degree 

options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt');
% options = optimoptions('fsolve','Jacobian','off');

for i=1:length(v:w)
    zz=v-1+i;
    fhandle=@(x)for_mainprog_3RPS(x,OS(:,zz));
    [xval(:,i),fval(:,i)] = fsolve(fhandle,x0,options);
    x0=xval(:,i);
    n1=xval(1,i);n3=xval(2,i);o3=xval(3,i);xc=xval(4,i);yc=xval(5,i);
    O1G(:,i)=[xc;yc;zc];   % w.r.t base coordinate system (bcs)
    O1P_bcs = Rot'*(OP-t);      % O1P_bcs = O1P described in bcs
    g = O1P_bcs-O1G(:,i);  % wrt bcs
    GP_bcs(:,i)=g/norm(g,2);    % wrt bcs
    GS=Rot'*OS(:,zz);         % Sun vector in bcs
    normal=(GS+GP_bcs(:,i))/norm((GS+GP_bcs(:,i)),2); % in bcs
%     GP_test(:,i) = (normal*norm((GS+GP_bcs(:,i)),2)-GS)/norm(((normal*norm((GS+GP_bcs(:,i)),2)-GS)),2);
    a1=normal(1,1);   a2=normal(2,1);       a3=normal(3,1);      % wrt base coordinate system
    n2=-yc/rp;
    o1=n2;
    o2=n1-(2*xc/rp);
    R(:,:,i) =     [n1 o1 a1 ;                     % rotation which takes mirror to base coordinate system
                    n2 o2 a2 ;
                    n3 o3 a3 ] ;
    Trans_matr(:,:,i)= [n1 o1 a1 xc;               % transformation which takes mirror to base coordinate system
                        n2 o2 a2 yc;
                        n3 o3 a3 zc;
                        0   0  0  1] ;
end  

TransmatrN=[Rot,t;[0 0 0 1]];           % transformation which takes base to global coordinate system

%% Just one mirror of dimension 2x2
s = 1.25;   % square receiver's half length
m = mirror_size/2;       % square mirror of half length 1. Actual mirror dimension is 2m x 2m.
PB1 = [s;s;0]; PB2 = [-s;s;0]; PB3 = [-s;-s;0]; PB4 = [s;-s;0];   % P is the centre of the receiver having coordinates [a;b;c]. 
                                                                  % B1, B2, B3, B4 are the corners of the square receiver
OB1 = Trans_Rec * [PB1;1]; OB2=Trans_Rec * [PB2;1]; OB3=Trans_Rec * [PB3;1]; OB4=Trans_Rec * [PB4;1];

Gm1= [m;m;0]; Gm2= [-m;m;0]; Gm3= [-m;-m;0]; Gm4= [m;-m;0];  % G is the centre of the mirror.
                                                             % m1, m2, m3 and m4 are the corners of the mirror

%%                                                             
figure()
mm=1; zz=0; 
options1 = optimoptions(@fmincon,'Algorithm','sqp');
% kg0 = [0.90;  -0.016;  -0.11;  0.70;  -0.10;  -0.30;  O1G(1,1);  O1G(2,1)];
kg0 = [0.90;  -0.016;  0.11;  0.014;  -0.1180;  0.90;  O1G(1,1);  O1G(2,1)];
% r0 = [0.1;-0.5;0.7];
for i=1:length(v:w)
    zz=v-1+i;
    f1 = @(kg)error_min(kg,a,b,c,OS(:,zz),zeta);
    f2 = @(kg)confuneq_error_min(kg);
%     [KG, fval] = fmincon(f1,kg0, [], [], [], [], [], [], f2,options1);
    [KG, fval] = fmincon(f1,kg0, [], [], [], [], [], [], f2);
    kg0=KG;
    KKG(:,i)=KG;
end

zz=0;
%%
for i=1:length(Trans_matr)
    zz=v-1+i;
    OG = t + Rot*O1G(:,i);  % in gcs
    Om1 = OG + R(:,:,i)*Gm1; Om2 = OG + R(:,:,i)*Gm2;  %% mirror corners in gcs
    Om3 = OG + R(:,:,i)*Gm3; Om4 = OG + R(:,:,i)*Gm4;
    
    % plot of receiver aperture
    l_Rec = plot3([OB1(1) OB2(1)], [OB1(2) OB2(2)], [OB1(3) OB2(3)],'c','LineWidth',3);
    hold on
    plot3([OB2(1) OB3(1)], [OB2(2) OB3(2)], [OB2(3) OB3(3)],'c','LineWidth',3)
    plot3([OB3(1) OB4(1)], [OB3(2) OB4(2)], [OB3(3) OB4(3)],'c','LineWidth',3)
    plot3([OB4(1) OB1(1)], [OB4(2) OB1(2)], [OB4(3) OB1(3)],'c','LineWidth',3)
    
    plot3([0 a],[0 b],[0,c],'b')     % Receiver tower
    plot3([0 a],[0 b],[0,c],'b*')     % Receiver tower center
    
%     % plot mirror
    plot3([Om1(1) Om2(1)], [Om1(2) Om2(2)], [Om1(3) Om2(3)],'g')
    plot3([Om2(1) Om3(1)], [Om2(2) Om3(2)], [Om2(3) Om3(3)],'g')
    plot3([Om3(1) Om4(1)], [Om3(2) Om4(2)], [Om3(3) Om4(3)],'g')
    plot3([Om4(1) Om1(1)], [Om4(2) Om1(2)], [Om4(3) Om1(3)],'g')
    plot3(OG(1),OG(2),OG(3),'g*')  % centroid of mirror
    
    % OHp = (Om1 + constant * unit reflected ray) should hit the receiver
    % aperture. OP is a point on the receiver surface where the sun ray
    % hitting the centre of the mirror supposed to fall. OHp - OP is a
    % vector in the plane of the receiver. The dot product of OHp and
    % normal to the receiver plane is zero. This is used to find the
    % constant
    
    GP = Rot*GP_bcs(:,i); % reflected ray in gcs
    k1=(cosd(zeta)*(a-Om1(1)) + sind(zeta)* (b-Om1(2)))/( GP(1)*cosd(zeta) + GP(2)*sind(zeta) );
    k2=(cosd(zeta)*(a-Om2(1)) + sind(zeta)* (b-Om2(2)))/( GP(1)*cosd(zeta) + GP(2)*sind(zeta) );
    k3=(cosd(zeta)*(a-Om3(1)) + sind(zeta)* (b-Om3(2)))/( GP(1)*cosd(zeta) + GP(2)*sind(zeta) );
    k4=(cosd(zeta)*(a-Om4(1)) + sind(zeta)* (b-Om4(2)))/( GP(1)*cosd(zeta) + GP(2)*sind(zeta) );
   
    % hitting points on the receiver
    hp1=Om1+k1*GP; % in gcs
    hp2=Om2+k2*GP; % in gcs
    hp3=Om3+k3*GP; % in gcs
    hp4=Om4+k4*GP; % in gcs
    centroid = (hp1+hp2+hp3+hp4)/4;
    plot3([0 centroid(1)], [0 centroid(2)], [0 centroid(3)], 'g*')
        
     % Plot of image on the receiver
    l_3RPS =   plot3([hp1(1) hp2(1)] ,[hp1(2) hp2(2)], [hp1(3) hp2(3)], 'g','LineWidth',1);
    plot3([hp2(1) hp3(1)] ,[hp2(2) hp3(2)], [hp2(3) hp3(3)], 'g','LineWidth',1)
    plot3([hp3(1) hp4(1)] ,[hp3(2) hp4(2)], [hp3(3) hp4(3)], 'g','LineWidth',1)
    plot3([hp4(1) hp1(1)] ,[hp4(2) hp1(2)], [hp4(3) hp1(3)], 'g','LineWidth',1)
    
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%     ERROR CORRECTION MODEL   %%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n1N = KKG(1,i); n3N = KKG(2,i); o3N = KKG(3,i);
    a1N = KKG(4,i);   a2N = KKG(5,i);       a3N = KKG(6,i);      % wrt base coordinate system
    xcN = KKG(7,i); ycN = KKG(8,i);
    n2N = -ycN/rp;
    o1N = n2N;
    o2N = n1N-(2*xcN/rp);
    O1GN = O1G(:,i);  % new G in bcs
    OGN = t+Rot*O1GN;   % This is the new point G in gcs
    RN  = [n1N o1N a1N ;               % Error corrected rotation matrix which 
           n2N o2N a2N ;               % takes mirror to base coordinate system
           n3N o3N a3N ];
    normalN = [a1N;a2N;a3N];  % in bcs
   
    Om1N = OGN + RN*Gm1; Om2N = OGN + RN*Gm2;  %% mirror corners in gcs
    Om3N = OGN + RN*Gm3; Om4N = OGN + RN*Gm4;
     % plot mirror
    plot3([Om1N(1) Om2N(1)], [Om1N(2) Om2N(2)], [Om1N(3) Om2N(3)],'r')
    plot3([Om2N(1) Om3N(1)], [Om2N(2) Om3N(2)], [Om2N(3) Om3N(3)],'r')
    plot3([Om3N(1) Om4N(1)], [Om3N(2) Om4N(2)], [Om3N(3) Om4N(3)],'r')
    plot3([Om4N(1) Om1N(1)], [Om4N(2) Om1N(2)], [Om4N(3) Om1N(3)],'r')
    plot3(OGN(1),OGN(2),OGN(3),'r*')  % centroid of mirror
%     
%     
    %% New reflected ray 
%     fhandle=@(r)Refl_ray(r,Rot'*OS(:,zz),normalN);
%     [GPN,fvald] = fsolve(fhandle,r0,options);
%     r0=GPN;
%     %%%%
    axy = cross(Rot'*OS(:,zz),[a1N;a2N;a3N]);  % IN bcs
    axs = axy/norm(axy,2);   % IN bcs
    angle1 = acos(dot(Rot'*OS(:,zz),[a1N;a2N;a3N]));  % IN bcs);  %%% RADIANS
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
%%%%%%%%%%%%%%%%%%%%
%%
    kon=5; kons=0.8;
    ref_ray = plot3([OGN(1) OGN(1)+kon*GPN(1)],[OGN(2) OGN(2)+kon*GPN(2)],[OGN(3) OGN(3)+kon*GPN(3)],'b-','LineWidth',2);
% %     ref_rayd = plot3([OGN(1) OGN(1)+kon*GPNd(1)],[OGN(2) OGN(2)+kon*GPNd(2)],[OGN(3) OGN(3)+kon*GPNd(3)],'g-','LineWidth',2);
% 
    sun = plot3([OGN(1) OGN(1)+kon*OS(1,zz)],[OGN(2) OGN(2)+kon*OS(2,zz)],[OGN(3) OGN(3)+kon*OS(3,zz)],'k-','LineWidth',2);
    nrml = plot3([OGN(1) OGN(1)+kon*normalN(1)],[OGN(2) OGN(2)+kon*normalN(2)],[OGN(3) OGN(3)+kon*normalN(3)],'r-','LineWidth',2);
    xaxis = plot3([OGN(1) OGN(1)+kons*RN(1,1)],[OGN(2) OGN(2)+kons*RN(2,1)],[OGN(3) OGN(3)+kons*RN(3,1)],'c-','LineWidth',2);
    yaxis = plot3([OGN(1) OGN(1)+kons*RN(1,2)],[OGN(2) OGN(2)+kons*RN(2,2)],[OGN(3) OGN(3)+kons*RN(3,2)],'m-','LineWidth',2);
    
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
    distN(i) = norm((centroidN-[a;b;c]),2);
    
    plot3([0 centroidN(1)], [0 centroidN(2)], [0 centroidN(3)], 'r*')
        
     % Plot of image on the receiver
    l_3RPSN =   plot3([hp1N(1) hp2N(1)] ,[hp1N(2) hp2N(2)], [hp1N(3) hp2N(3)], 'r','LineWidth',1);
    plot3([hp2N(1) hp3N(1)] ,[hp2N(2) hp3N(2)], [hp2N(3) hp3N(3)], 'r','LineWidth',1)
    plot3([hp3N(1) hp4N(1)] ,[hp3N(2) hp4N(2)], [hp3N(3) hp4N(3)], 'r','LineWidth',1)
    plot3([hp4N(1) hp1N(1)] ,[hp4N(2) hp1N(2)], [hp4N(3) hp1N(3)], 'r','LineWidth',1)
                       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    axis([tx-3 tx+3 ty-3 ty+3 0 tz+3]);   % to see the mirror 
%     axis([-1.5 1.5 -1.5 1.5 c-3 c+3]); %view(90,0)   % to see the Receiver
    view(90,0)
     axis square
%    axis equal
    xlabel('{\boldmath $ X, East$}')
    ylabel('{\boldmath $ Y, North$}')
    zlabel('{\boldmath $ Z, Zenith$}')
    grid on
    hold off
    M(mm)=getframe();
%     h1=legend([l_Rec l_AzEl l_TA],'Receiver','AzEl', 'TA',[300,67,3,4]);
%     set(h1,'Interpreter','latex');
%     M = getframe;
%     writeVideo(writerObj,M);
    mm=mm+1;
end
%%
dist = norm((OP-[a;b;c]));
