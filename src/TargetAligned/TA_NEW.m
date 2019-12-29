clear all
close all
home
set(0,'defaulttextinterpreter','latex')
% set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 15)
% set(0,'DefaultAxesFontWeight', 'Bold')

% Change default text fonts.
% set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 15)
% X axis-local east, Y axis - local north, Z axis - zenith (global coordinate system)

% Bangalore Latitude
lad=12;
lam=58;
las=13;
% Bangalore Longitude
lod=77;
lom=33;
los=37;

latitude=lad+(lam+las/60)/60;  % in degrees
longitude=lod+(lom+los/60)/60; % in degrees

year=2016;
month=3;
day=20;
d=date(year,month,day);

time_incre=120;  % seconds
LT=0*3600:time_incre:24*3600;

[elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude); % elevation and azimuth angles of sun
[xs ,ys, zs]=sunposition(elevation, azimuth);  % direction cosines of sun in the sky or sun vector
OS = [xs;ys;zs] ;  % Sun vector
a=0; b=0  ;  c=65;  % Receiver coordinates
OP = [a;b;c];   % Receiver vector from gcs origin, O
 
rad=100; % radius of the mirror from the origin
psi=30; % psi is the angle the heliostat centre makes with the local east(X) axis
% zG=0.2365;  % height of the heliostat
zG = 2;

% KJ = -0.0485; % distance between the mirror center, M, and the top platform center, G
KJ = -0.25;
O1M = [0;0;zG];   % vector from the base origin to the mirror origin described wrt bcs
OO1 = [rad*cosd(psi);rad*sind(psi);0];   % vector from the global origin to the base origin wrt gcs
 
lambda = atand(rad/(c-zG));  % angle the reflected ray makes with the vertical
gamma = 0; % the angle by which the x_b is inclined towards X (East)
Rot_gam=[ cosd(gamma) -sind(gamma) 0;
          sind(gamma)  cosd(gamma) 0;   % rotation of base CS about z axis
              0            0       1]; 
alpha = psi-gamma;          
R1 = [cosd(alpha) -sind(alpha) 0;
      sind(alpha)  cosd(alpha) 0;
          0             0      1]; 
R2 = [cosd(-lambda) 0  sind(-lambda);
           0        1       0       ;
      -sind(-lambda) 0 cosd(-lambda)] ;

R11 = R1*R2;
xx=R11(:,1); % direction perpendicular to reflected ray initially

% P is the receiver / target point
O1P_bcs=Rot_gam'*(OP-OO1);      % O1P_bcs = O1P described in bcs (base coordinate system)
g=O1P_bcs-O1M;
GP=g/norm(g,2);  % unit reflected ray wrt bcs

GP1 = (OP-OO1-Rot_gam*O1M)/norm((OP-OO1-Rot_gam*O1M),2); % unit reflected ray in global CS
v=8/time_incre*3600+1; % start time for the plots
w=9/time_incre*3600+1; % end time for the plots

jl = 200; % switch time
% ang = [zeros(1,jl),-10:-10:-10*(length(v:w)-jl)];
ang = -60;

for i=1:length(v:w)
    zz=v-1+i;
    %% for test
%     GS=OS(:,zz);         % Sun vector in global CS
%     GN(:,i)=(GS+GP1)/norm((GS+GP1),2);  % unit normal in global CS
    %%
   
    theta(i)=0.5*(acosd(dot(Rot_gam'*OS(:,zz),GP/norm(GP,2)))); % elevation angle
    X=dot(xx,Rot_gam'*OS(:,zz));
    Y=dot(cross(xx,Rot_gam'*OS(:,zz)),GP/norm(GP,2));
    rho(i)=atan2d(Y,X);   % spinning angle
    
   R3 = [cosd(rho(i)) -sind(rho(i)) 0;
         sind(rho(i))  cosd(rho(i)) 0;
             0             0        1];
   R4 = [cosd(theta(i))  0  sind(theta(i));
               0         1       0  ;
         -sind(theta(i)) 0 cosd(theta(i))] ;
   Rot_TA(:,:,i)  = R11 * R3 * R4;         % Rotation matrix at every time instant %%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    normal = [Rot_TA(1,3,i); Rot_TA(2,3,i); Rot_TA(3,3,i)];
    
    Skew_normal_matrix = [   0      -normal(3)    normal(2);
                         normal(3)     0        -normal(1);
                        -normal(2)   normal(1)       0   ];
    
    Rz(:,:,i) = expm((ang*pi/180)*Skew_normal_matrix);
    Ak = R1*R2;
    Rot_TA(:,:,i) = Rz(:,:,i)*Rot_TA(:,:,i);
%     Rot_TA(:,:,jl:end) = Rot_TAZ(:,:,jl:end);
%     Rot_TA(:,:,jl:end) = permute(Rot_TA(:,:,jl:end),[2,1,3]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%    Trans_matr(:,:,i) = [Rot_TA(:,:,i) O1M];  %% transformation which takes mirror to bcs
   O1G(:,i) = O1M+KJ*Rot_TA(:,3,i);   %% wrt bcs
end
%% Leg lengths
figure(1)
rb=0.50;  % circum radius of the base equilateral triangle
rp=0.50;  % circum radius of the platform

% ss=0.08;
% MM1=[ss;ss;0];  MM2=[-ss;ss;0];MM3=[-ss;-ss;0];  MM4=[ss;-ss;0];

O1R1=[rb;0;0];O1R2=[-.5*rb;sqrt(3)*rb/2;0];O1R3=[-.5*rb;-sqrt(3)*rb/2;0];
O1Rs=[O1R1 O1R2 O1R3];

GS1=[rp;0;0];GS2=[-.5*rp;sqrt(3)*rp/2;0];GS3=[-.5*rp;-sqrt(3)*rp/2;0];
GSs=[GS1 GS2 GS3];

for i=1:length(Rot_TA)
    
    O1S1 = Rot_TA(:,:,i)*GS1+O1G(:,i);   % w.r.t base coordinate system
    l1_dir(:,i)= O1S1-O1R1;  % direction of leg-length
%             l1_dir(:,i)= O1U1-O1R3;
    l1(:,i) = norm(l1_dir(:,i));   % leg length
        
    O1S2=Rot_TA(:,:,i)*GS2+O1G(:,i);
    l2_dir(:,i)=O1S2-O1R2;
    l2(:,i)=norm(l2_dir(:,i));
        
    O1S3=Rot_TA(:,:,i)*GS3+O1G(:,i);
    l3_dir(:,i)=O1S3-O1R3;
%             l3_dir(:,i)=O1U3-O1R1;
    l3(:,i)=norm(l3_dir(:,i));
end
plot(LT(v:w)/3600, l1*1000,'r','LineWidth', 2)
hold on
plot(LT(v:w)/3600, l2*1000,'g','LineWidth', 2)
plot(LT(v:w)/3600, l3*1000,'b','LineWidth', 2)
xlabel('Time, hrs')
ylabel(' leg length, mm')
% title('leg length Vs Time for T-A of 3-UPU')
legend({'leg1','leg2', 'leg3'},'interpreter','latex')

%% elevation and spinning angle plots
figure(2)
plot(LT(v:w)/3600 , theta, 'r','LineWidth', 2)
hold on
plot(LT(v:w)/3600 , rho, 'b','LineWidth', 2)
xlabel('Time, hrs')
ylabel('Angle, $^0$')
legend ('angle of incidence / elevation' , 'spinning angle')

%% For Animation
%  FigHandle = figure(3);
%   set(FigHandle, 'Position', [100, 100, 1049, 895]);
TransmatrN=[Rot_gam,OO1;[0 0 0 1]];           % transformation which takes base to global coordinate system
base=TransmatrN*[O1Rs;[1 1 1]];
% writerObj=VideoWriter('TA.mp4');
% writerObj.FrameRate = 12;
% open(writerObj);
m = 0.5;
MM1= [m;m;0]; MM2= [-m;m;0]; MM3= [-m;-m;0]; MM4= [m;-m;0];  % M is the centre of the mirror.

for i=1:length(Rot_TA)% 
%     for i=1:4
    xxx=base(1,:);                                   % wrt global coordinate system
    yyy=base(2,:);                                   % wrt global coordinate system
    zzz=base(3,:);                                   % wrt global coordinate system
    fill3(xxx,yyy,zzz,'r')
    hold on
    OM = OO1 + Rot_gam*O1M;
    Om1 = OO1 + Rot_gam*(O1M+Rot_TA(:,:,i)*MM1);
    Om2 = OO1 + Rot_gam*(O1M+Rot_TA(:,:,i)*MM2);
    Om3 = OO1 + Rot_gam*(O1M+Rot_TA(:,:,i)*MM3);
    Om4 = OO1 + Rot_gam*(O1M+Rot_TA(:,:,i)*MM4);
    
%     Mirror=[OM1 OM2 OM3 OM4];
%     xm=Mirror(1,:);
%     ym=Mirror(2,:);
%     zm=Mirror(3,:);
%     fill3(xm,ym,zm,'g')  % top platform or mirror
% Plot of mirror
%     plot3([Om1(1) Om2(1)] ,[Om1(2) Om2(2)], [Om1(3) Om2(3)], 'r','LineWidth',2)
%     plot3([Om2(1) Om3(1)] ,[Om2(2) Om3(2)], [Om2(3) Om3(3)], 'r','LineWidth',2)
%     plot3([Om3(1) Om4(1)] ,[Om3(2) Om4(2)], [Om3(3) Om4(3)], 'r','LineWidth',2)
%     plot3([Om4(1) Om1(1)] ,[Om4(2) Om1(2)], [Om4(3) Om1(3)], 'r','LineWidth',2)

    OS1=OO1+Rot_gam*(O1G(:,i)+Rot_TA(:,:,i)*GS1);
    OS2=OO1+Rot_gam*(O1G(:,i)+Rot_TA(:,:,i)*GS2);
    OS3=OO1+Rot_gam*(O1G(:,i)+Rot_TA(:,:,i)*GS3);
    top_triangle=[OS1 OS2 OS3];
    xt=top_triangle(1,:);
    yt=top_triangle(2,:);
    zt=top_triangle(3,:);
    fill3(xt,yt,zt,'c')  % top triangle
    
    plot3([xt(1) xxx(1)],[yt(1) yyy(1)],[zt(1) zzz(1)],'--k','LineWidth',2)
    plot3([xt(2) xxx(2)],[yt(2) yyy(2)],[zt(2) zzz(2)],'--k','LineWidth',2)
    plot3([xt(3) xxx(3)],[yt(3) yyy(3)],[zt(3) zzz(3)],'--k','LineWidth',2)
%     plot3([xt(1) xxx(3)],[yt(1) yyy(3)],[zt(1) zzz(3)],'--k','LineWidth',2.5)
%     plot3([xt(2) xxx(2)],[yt(2) yyy(2)],[zt(2) zzz(2)],'--k','LineWidth',2.5)
%     plot3([xt(3) xxx(1)],[yt(3) yyy(1)],[zt(3) zzz(1)],'--k','LineWidth',2.5)
    
    kons=1;kon=.5;zz=v-1+i;
%     normal=plot3([OM(1,1) OM(1,1)+kon*Rot_TA(1,3,i)],[OM(2,1) OM(2,1)+kon*Rot_TA(2,3,i)],[OM(3,1) OM(3,1)+kon*Rot_TA(3,3,i)],'r','LineWidth',2);
    
%     sunposition=plot3([OM(1,1) OM(1,1)+kons*xs(zz)],[OM(2,1) OM(2,1)+kons*ys(zz)],[OM(3,1) OM(3,1)+kons*zs(zz)],'k','LineWidth',2);   % position of the sun at every time_incre looking from plate coordinte system
%     plot3(OM(1,1)+kons*xs(v:zz),OM(2,1)+kons*ys(v:zz),OM(3,1)+kons*zs(v:zz),'k*');  % position of the sun at every time_incre looking from plate coordinte system
%     % %     coordintes of the receiver
%     plot3([OM(1,1) a],[OM(2,1) b],[OM(3,1) c],'b','LineWidth',2);
%     plot3([OM(1,1) a],[OM(2,1) b],[OM(3,1) c],'k*');
%     plot3([a a],[b b],[0 c],'b');
        kjs =.5;
%     plot3([OM(1,1) OM(1,1)+kjs*GP(1,1)],[OM(2,1) OM(2,1)+kjs*GP(2,1)],[OM(3,1) OM(3,1)+kjs*GP(3,1)],'b','LineWidth',2);
%     axis([tx-3 tx+3 ty-3 ty+3 0 1.5]); 
    
    %     xlim([10 30]); ylim([40 50]); zlim([0 20]);
    %     xlim([-50 120]); ylim([0 30]); zlim([0 50]);
    %     axis([-10 10 -5 5 60 70]); view(82,6)            % to see the
    %     Receiver
    
%     axis equal
    mg = 1;
    axis([OO1(1)-mg OO1(1)+mg OO1(2)-mg OO1(2)+mg 0 3.5]);  
%         axis([tx-1 tx+1 ty-1 ty+1 0 tz+3]);  view(24,20) % to see the mirror
%         view(-16,18) % to see the mirror
    axis square
    view(226,28);
    grid on
    xlabel('X, East','Rotation',-3)
    ylabel('Y, North ','Rotation',75)
    zlabel('Z, Zenith')
    hold off
    M(i)=getframe(gcf);
%         writeVideo(writerObj,M(i));
    %     pause(0.1)
end 
% close(writerObj)


  



