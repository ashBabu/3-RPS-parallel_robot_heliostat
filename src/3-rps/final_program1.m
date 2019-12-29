% this solves 5 equations in 5 unknowns
%% VERY IMP: for some months like june,  and some other months like december azimuth = 360-a(1). otherwise time vs az-El graph will not be continuous
%  As per the data sheet of chinese linear actuator( HF-TGE), the retracted length (when stroke is
%  0) is 1250mm. So the leg length should not be lower than this value.
%  From 1250mm, it can go to a max distance of 2250mm (1000 mm stroke). So
%  the maximum of each of the three legs should not exceed 2250 mm. To
%  obtain stow position, all the three legs to have equal length.

close all 
clear all
clc
home
set(0,'defaulttextinterpreter','latex')
% set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
% set(0,'DefaultAxesFontWeight', 'Bold')

% Change default text fonts.
% set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20)

% Bangalore Latitude
lad=12;
lam=58;
las=13;

% Bangalore Longitude
lod=77;
lom=33;
los=37;
% 
% latitude=lad+(lam+las/60)/60;  % in degrees
% longitude=lod+(lom+los/60)/60; % in degrees

%
% Rajasthan Coordinates
latitude = 27.0238;
longitude = 74.2179;

year=2016;
month=9;
day=22;
d=date(year,month,day);
% d=[79 141 266 356];
% d=79;
time_incre=60;
LT=0:time_incre:24*3600;

global OP rp zc t Rot
psi=30;
rad=100;
tx=rad*cosd(psi);ty=rad*sind(psi);tz=0; t=[tx;ty;tz];
a=0;b=0;c=65;
% rb=0.36;
rb = 0.50;
rp=0.50;
[elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude);
OP=[a;b;c];
[xs ,ys, zs]=sunposition(elevation, azimuth);  % xs,ys,zs are coordinates of the sun
OS=[xs;ys;zs];
zc=2;
O1R1=[rb;0;0];O1R2=[-.5*rb;sqrt(3)*rb/2;0];O1R3=[-.5*rb;-sqrt(3)*rb/2;0];
O1Rs=[O1R1 O1R2 O1R3];

GS1=[rp;0;0];GS2=[-.5*rp;sqrt(3)*rp/2;0];GS3=[-.5*rp;-sqrt(3)*rp/2;0];
GSs=[GS1 GS2 GS3];
%%
gamma=0;
Rot=[ cosd(gamma) -sind(gamma) 0;
      sind(gamma)  cosd(gamma) 0;
          0            0       1]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x0 = [.7;-.6;0.1;-0.7;-0.01]; % works for psi=0,45,90,135,180,225,270,315 degrees
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x0 = [0.2;0.61;0.5;-0.01;0.015];  % works only for psi=45 degree 
        
% Some other tested values
% x0=[0.25; 0.6; 0.4; 0.35; 0.56];
% x0 = [1.7;.6;-1;.7;.1;.5;0.08;1]; 
x0 = [0.1761;-0.5876;-0.7897;0.9579;0.2871];

options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt');
% options = optimoptions('fsolve','Display','iter','TolX',1e-9);
% options = optimoptions('fsolve','Jacobian','off')
% base1=base_platform(alpha,rb);           % wrt base coordinate system
% top=top_platform(beta,rp);              % wrt mirror coordinate system
v=7/time_incre*3600+1;
w=18/time_incre*3600+1;

%% Plot of azimuth and elevation angles
figure(1)
plot(LT(v:w)/3600,elevation((v: w)),'r','LineWidth',2)
hold on
plot(LT(v:w)/3600,azimuth((v: w)),'b','LineWidth',2)
xlabel('Time, hrs')
ylabel('Angle, Degrees')
% title('Variation of the leg lengths w.r.t. time')
lg = legend({' Elevation',' Azimuth'},'interpreter','latex');
%set(lg,'defaulttextinterpreter','latex')
grid on
%%
for i=1:length(v:w)
    zz=v-1+i;
    fhandle=@(x)for_mainprog_3RPS(x,OS(:,zz));
    [xval(:,i),fval(:,i)] = fsolve(fhandle,x0,options);
    x0=xval(:,i);
    n1=xval(1,i);n3=xval(2,i);o3=xval(3,i);xc=xval(4,i);yc=xval(5,i);
    B0t(:,i)=[xc;yc;zc];                        
    O1G=[xc;yc;zc];   % w.r.t base coordinate system (bcs)
    O1P_bcs=Rot'*(OP-t);      % O1P_bcs = O1P described in bcs
    g=O1P_bcs-O1G;
    GP(:,i)=g/norm(g,2);
    GS=Rot'*OS(:,zz);         % Sun vector in bcs
    normal=(GS+GP(:,i))/norm((GS+GP(:,i)),2);
    a1=normal(1,1);   a2=normal(2,1);       a3=normal(3,1);      % wrt base coordinate system
    n2=-yc/rp;
    o1=n2;
    o2=n1-(2*xc/rp);
    R(:,:,i) =     [n1 o1 a1 ;                               % rotation which takes mirror to base coordinate system
                    n2 o2 a2 ;
                    n3 o3 a3 ] ;
%     Trans_matr(:,:,i)=[R(:,:,i), [xc;yc;zc]; [0  0   0  1  ] ]; 
    Trans_matr(:,:,i)= [n1 o1 a1 xc;                               % transformation which takes mirror to base coordinate system
                        n2 o2 a2 yc;
                        n3 o3 a3 zc;
                        0   0  0  1] ;
end      

%% Variation of the centre of the platform
figure(2)
plot(B0t(1,:),B0t(2,:),'r','LineWidth',2)
grid on
hold on
xlabel('X, m')
ylabel('Y, m')
xlim([-0.16 0.1])
ylim([-0.05 0.15])

aw1 = 9/time_incre*3600+1; aw2 = 12/time_incre*3600+1; aw3 = 15/time_incre*3600+1; aw4 = 17/time_incre*3600+1;
txt1 = (['T =  ' num2str((aw1-1)*time_incre/3600) ':00' ]);
txt2 = (['T =  ' num2str((aw2-1)*time_incre/3600) ':00' ]);
txt3 = (['T =  ' num2str((aw3-1)*time_incre/3600) ':00' ]);
txt4 = (['T =  ' num2str((aw4-1)*time_incre/3600) ':00' ]);
% legend([p1 p2 p3],txt1,txt2,txt3,'Location','northwest', 'Interpreter', 'latex');
set(gca,'TickLabelInterpreter', 'latex');
% txt1 = (['$\psi$ = ' num2str(psi) '$^0$']);
plot(B0t(1,aw1-v+1),B0t(2,aw1-v+1),'*', 'MarkerSize', 20)
text(B0t(1,aw1-v+1)-0.08,B0t(2,aw1-v+1),txt1)
plot(B0t(1,aw2-v+1),B0t(2,aw2-v+1),'*', 'MarkerSize', 20)
text(B0t(1,aw2-v+1)-0.085,B0t(2,aw2-v+1),txt2)
plot(B0t(1,aw3-v+1),B0t(2,aw3-v+1),'*', 'MarkerSize', 20)
text((B0t(1,aw3-v+1))-0.018,B0t(2,aw3-v+1)+0.015,txt3)
plot(B0t(1,aw4-v+1),B0t(2,aw4-v+1),'*', 'MarkerSize', 20)
text((B0t(1,aw4-v+1))-0.05,B0t(2,aw4-v+1)-0.02,txt4)
% zlabel('{\boldmath $ Z, m$}')
%%
    
TransmatrN=[Rot,t;[0 0 0 1]];           % transformation which takes base to global coordinate system
base=TransmatrN*[O1Rs;[1 1 1]];

%% 
for i=1:length(Trans_matr)
    if i==1
        ad=1;
    else
        ad=i-1;
    end
    B0p1(:,:,i)=Trans_matr(:,:,i)*[GS1;1];   % w.r.t base coordinate system
    l1_dir(:,i)=B0p1(:,:,i)-[O1R1;1];
    l1(i,:)=norm(l1_dir(:,i));
%     vel1(:,i)=(l1(:,i) - l1(:,ad));%/time_incre;
    theta1(i,:)=180-acosd(dot(O1R1,l1_dir(1:3,i))/(norm(O1R1,2)*l1(i,:)));
    
    B0p2(:,:,i)=(Trans_matr(:,:,i)*[GS2;1]);
    l2_dir(:,i)=(B0p2(:,:,i)-[O1R2;1]);
    l2(i,:)=norm(l2_dir(:,i));
%     vel2(:,i)=(l2(:,i) - l2(:,ad));%/time_incre;
    theta2(i,:)=180-acosd(dot(O1R2,l2_dir(1:3,i))/(norm(O1R2,2)*norm(l2_dir(:,i),2)));
    
    B0p3(:,:,i)=(Trans_matr(:,:,i)*[GS3;1]);
    l3_dir(:,i)=(B0p3(:,:,i)-[O1R3;1]);
    l3(i,:)=norm(l3_dir(:,i));
%     vel3(:,i)=(l3(i,:) - l3(ad,:));%/time_incre;
    theta3(i,:)=180-acosd(dot(O1R3,l3_dir(1:3,i))/(norm(O1R3,2)*norm(l3_dir(:,i),2)));
end
%% plot of leg length
figure (3)
plot(LT(v:w)/3600,1000*l1(1:length(v: w)),'r--o','LineWidth',2)
hold on
plot(LT(v:w)/3600,1000*l2(1:length(v: w)),'g+','LineWidth',2)
plot(LT(v:w)/3600,1000*l3(1:length(v: w)),'b','LineWidth',2)

xlabel('Time, hrs')
ylabel('Leg length, mm')
% title('leg length variation wrt time')
legend({'leg1','leg2', 'leg3'},'interpreter','latex')
grid on
hold on

% stroke_l1=max(l1)-min(l1)
% stroke_l2=max(l2)-min(l2)
% stroke_l3=max(l3)-min(l3
% Min_vel_leg1=mean(abs(vel1(2:end)))
% Min_vel_leg2=mean(abs(vel2(2:end)))
% Min_vel_leg3=mean(abs(vel3(2:end)))
%% theta plot
 figure(4)
 plot(LT(v:w)/3600,90-theta1(1:length(v: w)),'g','LineWidth',2)
 hold on
 plot(LT(v:w)/3600,90-theta2(1:length(v: w)),'r','LineWidth',2)
 plot(LT(v:w)/3600,90-theta3(1:length(v: w)),'b','LineWidth',2)
 
 xlabel('Time, hrs')
ylabel('Angle, Degrees')
% % title('Variation of the thetas from vertical w.r.t. time')
legend({'$\theta_1$','$\theta_2$','$\theta_3$'},'interpreter','latex')
grid on   

%% Animation
figure(5)
%   FigHandle = figure(5);
%   set(FigHandle, 'Position', [100, 100, 1049, 895]);
% writerObj=VideoWriter('3RPS_field.avi');
% writerObj.FrameRate = 6;
% open(writerObj);
for i=1:length(v:w)   % from 7 a.m to 6 p.m (18:00))
    % for i=1:4
    xxx=base(1,:);                                   % wrt global coordinate system
    yyy=base(2,:);                                   % wrt global coordinate system
    zzz=base(3,:);                                   % wrt global coordinate system
    fill3(xxx,yyy,zzz,'r')
    hold on
    OG=Rot*B0t(:,i)+t;
    OS1=OG+Rot*R(:,:,i)*GS1;
    OS2=OG+Rot*R(:,:,i)*GS2;
    OS3=OG+Rot*R(:,:,i)*GS3;
    top=[OS1 OS2 OS3];
    x=top(1,:);
    y=top(2,:);
    z=top(3,:);
    fill3(x,y,z,'g')  % top platform or mirror
    TRANSMATR=Rot*R(:,:,i);
    %     %% #1 change the leg lengths also if the below plot is changed
    plot3([x(1) xxx(1)],[y(1) yyy(1)],[z(1) zzz(1)],'--k','LineWidth',2.5)
    plot3([x(2) xxx(2)],[y(2) yyy(2)],[z(2) zzz(2)],'--k','LineWidth',2.5)
    plot3([x(3) xxx(3)],[y(3) yyy(3)],[z(3) zzz(3)],'--k','LineWidth',2.5)
    %    %% changing wrt global coordinate system
    
    kons=20;kon=20;zz=v-1+i;
    normal = plot3([OG(1,1) OG(1,1)+kon*TRANSMATR(1,3)],[OG(2,1) OG(2,1)+kon*TRANSMATR(2,3)],[OG(3,1) OG(3,1)+kon*TRANSMATR(3,3)],'r:','LineWidth',2);
    %     proj=plot3([B0t(1,i) B0t(1,i)+kon*Trans_matr(1,3,i)],[B0t(2,i) B0t(2,i)+kon*Trans_matr(2,3,i)],[0 0],'r','LineWidth',2);
    
    sunposition=plot3([OG(1,1) OG(1,1)+kons*xs(zz)],[OG(2,1) OG(2,1)+kons*ys(zz)],[OG(3,1) OG(3,1)+kons*zs(zz)],'--k','LineWidth',2);   % position of the sun at every time_incre looking from plate coordinte system
    plot3(OG(1,1)+kons*xs(v:zz),OG(2,1)+kons*ys(v:zz),OG(3,1)+kons*zs(v:zz),'k*');  % position of the sun at every time_incre looking from plate coordinte system
    % %     coordintes of the receiver
    GP = (OP-OG)/norm((OP-OG),2);  % Unit reflected ray
    ref_ray = plot3([OG(1,1) OG(1,1)+kon*GP(1)],[OG(2,1) OG(2,1)+kon*GP(2)],[OG(3,1)  OG(3,1)+kon*GP(3)],'b','LineWidth',2);

    plot3([OG(1,1) a],[OG(2,1) b],[OG(3,1) c],'b','LineWidth',2);
    plot3([OG(1,1) a],[OG(2,1) b],[OG(3,1) c],'k*');
    plot3([a a],[b b],[0 c],'b');
        axis equal
    %     xlim([10 30]); ylim([40 50]); zlim([0 20]);
    %     xlim([-50 120]); ylim([0 30]); zlim([0 50]);
        axis([tx-1 tx+1 ty-1 ty+1 0 tz+3]);  view(24,20) % to see the mirror
%         axis([-10 10 -5 5 60 70]); view(82,6)            % to see the
    %     Receiver
%     axis square

%     view(-26,34);
%     view(az_n(i),el_n(i))
    grid on
%     xlabel('X, East')
%     ylabel('Y, North')
%      xlabel(' X, East','Rotation',28,'Position',[1167.5 1416 -1048.7])
%     ylabel('Y, North','Rotation',-35,'Position',[1087.5 1473.8 -1053.1])
     xlabel(' X, East','Rotation',-5)
    ylabel('Y, North','Rotation',60)
    zlabel(' Z, Zenith')
    legend([normal,sunposition,ref_ray],'Normal','Incident ray','Reflected ray')
    hold off
     M(i)=getframe(gcf);
%         writeVideo(writerObj,M(i));
%         pause(0.1)
end
a
% close(writerObj)

% writerObj=VideoWriter('solar1.avi');
% writerObj.FrameRate = 9;
% open(writerObj);
%%
% for i=1:length(v:w)
%       xn=Trans_matr(1,3,i)  ;      yn=Trans_matr(2,3,i)  ;     zn=Trans_matr(3,3,i)  ;   prj=[xn;yn;0];    y_axis=[0;1;0];
%       ele_n(i)=atan2d(zn,norm(prj,2));
%       flag = atan2d(yn,xn);
%       if (flag<0 && flag>=-360)
%           az_n(i)=rem(360+flag,360);
%       elseif (flag<-360)
%           az_n(i)=mod(360+flag,360);
%       else
%           az_n(i)=flag;
%       end
% end
% figure(5)
% plot(LT(v:w)/3600,ele_n,'r*','LineWidth',4)
% hold on
% plot(LT(v:w)/3600,az_n,'b*','LineWidth',4)
% grid on
% xlabel('Time  \it{t / hr}')
% ylabel('Angles subtended by heliostat normal  \it{\theta / (^0)}' )
% legend('Elevation','Azimuth')
% 
% 


% %% velocity plot
%  figure (6)
% plot(LT(v:w)/3600,vel1(1:length(v: w)),'r*','LineWidth',4)
% hold on
% plot(LT(v:w)/3600,vel2(1:length(v: w)),'b*','LineWidth',4)
% plot(LT(v:w)/3600,vel3(1:length(v: w)),'g*','LineWidth',4)
% % xlabel('Time  \it{t / hr}')
% % ylabel('Leg lengths  \it{ l / m} ')
% % title('Variation of the leg lengths w.r.t. time')
% lgnd=legend('velocity l_1','velocity l_2','velocity l_3');
% grid on
% 
% %%
% [A1, A2, A3,B1,B2,B3]=XYZ_Euler_invert(R,theta1,theta2,theta3);
% [A11, A21, A31,B11,B21,B31]=ZYZ_Euler_invert(R,theta1,theta2,theta3);   
[A10, A20, A30,B10,B20,B30]=ZYX_Euler_invert(R,theta1,theta2,theta3);
% %% Values obtained from ADAMS
% time = [0,504,1010,1510,2020,2520,3020,3530,4030,4540,5040,5540,6050,6550,7060,7560,8060,8570,9070,9580,10100,10600,11100,11600,12100,12600,13100,13600,14100,14600,15100,15600,16100,16600,17100,17600,18100,18600,19200,19700,20200,20700,21200,21700,22200,22700,23200,23700,24200,24700,25200];
% rb36_rp5_SpherRot_Leg1_X = [0.262000000000000,-0.751000000000000,-1.79000000000000,-2.84000000000000,-3.93000000000000,-5.03000000000000,-6.16000000000000,-7.30000000000000,-8.46000000000000,-9.63000000000000,-10.8000000000000,-12,-13.3000000000000,-14.5000000000000,-15.7000000000000,-17,-18.3000000000000,-19.6000000000000,-20.9000000000000,-22.2000000000000,-23.5000000000000,-24.8000000000000,-26.1000000000000,-27.5000000000000,-28.8000000000000,-30.2000000000000,-31.5000000000000,-32.9000000000000,-34.2000000000000,-35.6000000000000,-37,-38.3000000000000,-39.7000000000000,-41.1000000000000,-42.4000000000000,-43.8000000000000];
% rb36_rp5_SpherRot_Leg1_Y = [5.42000000000000,4.39000000000000,3.41000000000000,2.46000000000000,1.55000000000000,0.660000000000000,-0.192000000000000,-1.02000000000000,-1.82000000000000,-2.59000000000000,-3.33000000000000,-4.05000000000000,-4.74000000000000,-5.42000000000000,-6.07000000000000,-6.69000000000000,-7.29000000000000,-7.86000000000000,-8.42000000000000,-8.96000000000000,-9.47000000000000,-9.96000000000000,-10.4000000000000,-10.9000000000000,-11.3000000000000,-11.7000000000000,-12.1000000000000,-12.4000000000000,-12.8000000000000,-13.1000000000000,-13.4000000000000,-13.6000000000000,-13.9000000000000,-14.1000000000000,-14.3000000000000,-14.4000000000000];
% 
% rb36_rp5_SpherRot_Leg2_X = [6.90000000000000,6.52000000000000,6.20000000000000,5.91000000000000,5.67000000000000,5.46000000000000,5.29000000000000,5.15000000000000,5.05000000000000,4.97000000000000,4.92000000000000,4.90000000000000,4.91000000000000,4.93000000000000,4.99000000000000,5.08000000000000,5.18000000000000,5.32000000000000,5.47000000000000,5.64000000000000,5.85000000000000,6.09000000000000,6.34000000000000,6.62000000000000,6.94000000000000,7.27000000000000,7.65000000000000,8.05000000000000,8.48000000000000,8.97000000000000,9.47000000000000,10,10.6000000000000,11.3000000000000,12,12.7000000000000];
% rb36_rp5_SpherRot_Leg2_Y = [-6.91000000000000,-5.57000000000000,-4.22000000000000,-2.86000000000000,-1.48000000000000,-0.0908000000000000,1.31000000000000,2.73000000000000,4.16000000000000,5.59000000000000,7.04000000000000,8.49000000000000,9.95000000000000,11.4000000000000,12.9000000000000,14.4000000000000,15.9000000000000,17.4000000000000,18.9000000000000,20.4000000000000,21.9000000000000,23.4000000000000,24.9000000000000,26.5000000000000,28,29.5000000000000,31,32.5000000000000,34.1000000000000,35.6000000000000,37.1000000000000,38.6000000000000,40.1000000000000,41.5000000000000,43,44.4000000000000];

% rb36_rp5_SpherRot_Leg3_X = [-7.16000000000000,-5.78000000000000,-4.43000000000000,-3.08000000000000,-1.75000000000000,-0.435000000000000,0.872000000000000,2.18000000000000,3.46000000000000,4.73000000000000,6.01000000000000,7.26000000000000,8.51000000000000,9.76000000000000,11,12.2000000000000,13.5000000000000,14.7000000000000,15.9000000000000,17.1000000000000,18.3000000000000,19.6000000000000,20.8000000000000,22,23.2000000000000,24.4000000000000,25.6000000000000,26.7000000000000,27.9000000000000,29.1000000000000,30.3000000000000,31.5000000000000,32.7000000000000,33.8000000000000,35,36.2000000000000];
% rb36_rp5_SpherRot_Leg3_Y = [-6.46000000000000,-6.83000000000000,-7.23000000000000,-7.67000000000000,-8.14000000000000,-8.64000000000000,-9.17000000000000,-9.73000000000000,-10.3000000000000,-10.9000000000000,-11.5000000000000,-12.2000000000000,-12.8000000000000,-13.5000000000000,-14.2000000000000,-14.9000000000000,-15.6000000000000,-16.3000000000000,-17,-17.7000000000000,-18.5000000000000,-19.2000000000000,-19.9000000000000,-20.7000000000000,-21.4000000000000,-22.2000000000000,-22.9000000000000,-23.7000000000000,-24.4000000000000,-25.2000000000000,-25.9000000000000,-26.7000000000000,-27.4000000000000,-28.2000000000000,-28.9000000000000,-29.6000000000000];


% %% spherical joint angle plot
% figure (7) % Rotaion of 1st spherical joint on the xz plane
% hFig1 = figure(7);
% set(hFig1, 'Position', [0 0 1800 1800])
% % mtit(hFig1,'XYZ Euler Inversion angles')
% subplot(2,3,1); 
% plot(LT(v:w)/3600,A1(1,:),'r*','LineWidth',2)   % X rot leg1
% hold on
% plot(LT(v:w)/3600,A1(2,:),'b*','LineWidth',2)   % Y rot leg1
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_X ,'r','LineWidth',1 )  % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_Y ,'b','LineWidth',1 )  % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 1')
% legend('Rot abt X',' Rot abt Y')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (8) % Rotaion of 2nd spherical joint on the y = -sqrt(3) plane
% subplot(2,3,2); 
% plot(LT(v:w)/3600,A2(1,:),'r*','LineWidth',2)       % X rot leg2
% hold on
% plot(LT(v:w)/3600,A2(2,:),'b*','LineWidth',2)       % Y rot leg2
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 2')
% legend('Rot abt X',' Rot abt Y')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (9) % Rotaion of 3rd spherical joint on the y = sqrt(3) plane
%  subplot(2,3,3);
% plot(LT(v:w)/3600,A3(1,:),'r*','LineWidth',2)       % X rot leg3
% hold on
% plot(LT(v:w)/3600,A3(2,:),'b*','LineWidth',2)       % Y rot leg3
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% legend('Rot abt X',' Rot abt Y')
% title('Variation of the spherical joint angles w.r.t. time for leg 3')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% % figure (10) % Rotaion of 1st spherical joint on the xz plane
% subplot(2,3,4); 
% plot(LT(v:w)/3600,B1(1,:),'r*','LineWidth',2)   % X rot leg1
% hold on
% plot(LT(v:w)/3600,B1(2,:),'b*','LineWidth',2)   % Y rot leg1
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_X ,'r','LineWidth',1 )  % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_Y ,'b','LineWidth',1 )  % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 1')
% legend('Rot abt X',' Rot abt Y')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (11) % Rotaion of 2nd spherical joint on the y = -sqrt(3) plane
% subplot(2,3,5); 
% plot(LT(v:w)/3600,B2(1,:),'r*','LineWidth',2)       % X rot leg2
% hold on
% plot(LT(v:w)/3600,B2(2,:),'b*','LineWidth',2)       % Y rot leg2
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 2')
% legend('Rot abt X',' Rot abt Y')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (12) % Rotaion of 3rd spherical joint on the y = sqrt(3) plane
% subplot(2,3,6); 
% plot(LT(v:w)/3600,B3(1,:),'r*','LineWidth',2)       % X rot leg3
% hold on
% plot(LT(v:w)/3600,B3(2,:),'b*','LineWidth',2)       % Y rot leg3
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% legend('Rot abt X',' Rot abt Y')
% title('Variation of the spherical joint angles w.r.t. time for leg 3')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% % mtit(hFig1,'XYZ Euler Inversion angles')
% suptitle('XYZ Euler Inversion angles')
% 
% 
% %% ZYZ Euler Inversions
% 
% figure (8) % Rotaion of 1st spherical joint on the xz plane
% hFig = figure(8);
% set(hFig, 'Position', [0 0 1800 1800])
% 
% subplot(2,3,1); 
% plot(LT(v:w)/3600,A11(1,:),'r*','LineWidth',2)   % Z rot leg1
% hold on
% plot(LT(v:w)/3600,A11(2,:),'b*','LineWidth',2)   % Y rot leg1
% plot(LT(v:w)/3600,A11(3,:),'g*','LineWidth',2)   % Z* rot leg1
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_X ,'r','LineWidth',1 )  % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_Y ,'b','LineWidth',1 )  % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 1')
% legend('Rot abt Z',' Rot abt Y','Rot abt Z*')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (8) % Rotaion of 2nd spherical joint on the y = -sqrt(3) plane
% subplot(2,3,2); 
% plot(LT(v:w)/3600,A21(1,:),'r*','LineWidth',2)       % Z rot leg2
% hold on
% plot(LT(v:w)/3600,A21(2,:),'b*','LineWidth',2)       % Y rot leg2
% plot(LT(v:w)/3600,A21(3,:),'g*','LineWidth',2)       % Z* rot leg2
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 2')
% legend('Rot abt Z',' Rot abt Y', 'Rot abt Z*')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (9) % Rotaion of 3rd spherical joint on the y = sqrt(3) plane
%  subplot(2,3,3);
% plot(LT(v:w)/3600,A31(1,:),'r*','LineWidth',2)       % Z rot leg3
% hold on
% plot(LT(v:w)/3600,A31(2,:),'b*','LineWidth',2)       % Y rot leg3
% plot(LT(v:w)/3600,A31(3,:),'g*','LineWidth',2)       % Z* rot leg3
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% legend('Rot abt Z',' Rot abt Y','Rot abt Z*')
% title('Variation of the spherical joint angles w.r.t. time for leg 3')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% % figure (10) % Rotaion of 1st spherical joint on the xz plane
% subplot(2,3,4); 
% plot(LT(v:w)/3600,B11(1,:),'r*','LineWidth',2)   % Z rot leg1
% hold on
% plot(LT(v:w)/3600,B11(2,:),'b*','LineWidth',2)   % Y rot leg1
% plot(LT(v:w)/3600,B11(3,:),'g*','LineWidth',2)   % Z* rot leg1
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_X ,'r','LineWidth',1 )  % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg1_Y ,'b','LineWidth',1 )  % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 1')
% legend('Rot abt Z',' Rot abt Y','Rot abt Z*')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (11) % Rotaion of 2nd spherical joint on the y = -sqrt(3) plane
% subplot(2,3,5); 
% plot(LT(v:w)/3600,B21(1,:),'r*','LineWidth',2)       % Z rot leg2
% hold on
% plot(LT(v:w)/3600,B21(2,:),'b*','LineWidth',2)       % Y rot leg2
% plot(LT(v:w)/3600,B21(3,:),'g*','LineWidth',2)       % Z* rot leg2
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg2_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 2')
% legend('Rot abt Z',' Rot abt Y','Rot abt Z*')
% % legend('MATLAB Rot abt X',' MATLAB Rot abt Y', 'ADAMS Rot abt X',' ADAMS Rot abt Y')
% % grid on
% 
% %  figure (12) % Rotaion of 3rd spherical joint on the y = sqrt(3) plane
% subplot(2,3,6); 
% plot(LT(v:w)/3600,B31(1,:),'r*','LineWidth',2)       % Z rot leg3
% hold on
% plot(LT(v:w)/3600,B31(2,:),'b*','LineWidth',2)       % Y rot leg3
% plot(LT(v:w)/3600,B31(3,:),'g*','LineWidth',2)       % Z* rot leg3
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_X ,'r','LineWidth',1 )   % Obtained from adams
% % plot(LT(v:w)/3600 , rb36_rp5_SpherRot_Leg3_Y ,'b','LineWidth',1 )   % Obtained from adams
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% legend('Rot abt Z',' Rot abt Y','Rot abt Z*')
% title('Variation of the spherical joint angles w.r.t. time for leg 3')
% suptitle('ZYZ Euler Inversion angles')

%% ZYX Euler Inverstion
% figure (9) % Rotaion of 1st spherical joint on the xz plane
% hFig = figure(9);
% set(hFig, 'Position', [0 0 1800 1800])
% 
% subplot(2,3,1); 
% plot(LT(v:w)/3600,A10(1,:),'r*','LineWidth',2)   % Z rot leg1
% hold on
% plot(LT(v:w)/3600,A10(2,:),'b*','LineWidth',2)   % Y rot leg1
% plot(LT(v:w)/3600,A10(3,:),'g*','LineWidth',2)   % X rot leg1
% 
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 1')
% legend('Rot abt Z',' Rot abt Y','Rot abt X')
% 
% subplot(2,3,2); 
% plot(LT(v:w)/3600,A20(1,:),'r*','LineWidth',2)       % Z rot leg2
% hold on
% plot(LT(v:w)/3600,A20(2,:),'b*','LineWidth',2)       % Y rot leg2
% plot(LT(v:w)/3600,A20(3,:),'g*','LineWidth',2)       % X rot leg2
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 2')
% legend('Rot abt Z',' Rot abt Y', 'Rot abt X')
% 
% subplot(2,3,3);
% plot(LT(v:w)/3600,A30(1,:),'r*','LineWidth',2)       % Z rot leg3
% hold on
% plot(LT(v:w)/3600,A30(2,:),'b*','LineWidth',2)       % Y rot leg3
% plot(LT(v:w)/3600,A30(3,:),'g*','LineWidth',2)       % X rot leg3
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% legend('Rot abt Z',' Rot abt Y','Rot abt X')
% title('Variation of the spherical joint angles w.r.t. time for leg 3')
% 
% subplot(2,3,4); 
% plot(LT(v:w)/3600,B10(1,:),'r*','LineWidth',2)   % Z rot leg1
% hold on
% plot(LT(v:w)/3600,B10(2,:),'b*','LineWidth',2)   % Y rot leg1
% plot(LT(v:w)/3600,B10(3,:),'g*','LineWidth',2)   % X rot leg1
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 1')
% legend('Rot abt Z',' Rot abt Y','Rot abt X')
% 
% subplot(2,3,5); 
% plot(LT(v:w)/3600,B20(1,:),'r*','LineWidth',2)       % Z rot leg2
% hold on
% plot(LT(v:w)/3600,B20(2,:),'b*','LineWidth',2)       % Y rot leg2
% plot(LT(v:w)/3600,B20(3,:),'g*','LineWidth',2)       % X rot leg2
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% title('Variation of the spherical joint angles w.r.t. time for leg 2')
% legend('Rot abt Z',' Rot abt Y','Rot abt X')
% 
% subplot(2,3,6); 
% plot(LT(v:w)/3600,B30(1,:),'r*','LineWidth',2)       % Z rot leg3
% hold on
% plot(LT(v:w)/3600,B30(2,:),'b*','LineWidth',2)       % Y rot leg3
% plot(LT(v:w)/3600,B30(3,:),'g*','LineWidth',2)       % X rot leg3
% xlabel('Time  \it{t / hr}')
% ylabel('Angle  \it{ \theta / ^0} ')
% legend('Rot abt Z',' Rot abt Y','Rot abt X')
% title('Variation of the spherical joint angles w.r.t. time for leg 3')
% str = sprintf(' ZYX Euler Inversion angles: radius = %4.1f and \\psi = %2.1f',rad,psi);
% % str = text(0.5,1,' ZYX Euler Inversion angles: radius = %4.1f and \theta psi = %2.1f',rad,psi);
% suptitle(str);
% suptitle('ZYX Euler Inversion angles', )


%% To make a txt file to give as input spline to ADAMS
% BL_Cylin_Radius= 0.2; SL_Cylin_Leng=2.3; Sphere_Dia=BL_Cylin_Radius-0.1;
% leg1 = fopen('first_leg.txt','w');
% g=0;
% formatSpec = ' %4.1f  \t %4.4f \n';
% for i=1:length(v:w)
%     fprintf(leg1,formatSpec,g,l1(i)-SL_Cylin_Leng-Sphere_Dia)
%     g=g+time_incre;
% end
% table1='first_leg.txt';
% fclose(leg1);
% 
% leg2 = fopen('second_leg.txt','w');
% g=0;
% formatSpec = ' %4.1f  \t %4.4f \n';
% for i=1:length(v:w)
%     fprintf(leg2,formatSpec,g,l2(i)-SL_Cylin_Leng-Sphere_Dia)
%     g=g+time_incre;
% end
% table2='second_leg.txt';
% fclose(leg2);
% 
% leg3 = fopen('third_leg.txt','w');
% g=0;
% formatSpec = ' %4.1f  \t %4.4f \n';
% for i=1:length(v:w)
%     fprintf(leg3,formatSpec,g,l3(i)-SL_Cylin_Leng-Sphere_Dia)
%     g=g+time_incre;
% end
% table3='third_leg.txt';
% fclose(leg3);
%%  The below figures are repeatitions of the big figure(subplots).
% % This was done for the paper for solarpaces. 
% The solarpaces final paper Rot abt X and Rot abt Y graphs are wrong

figure(15)
% plot(LT(v:w)/3600,A10(1,:),'r*','LineWidth',2)   % Z rot leg1
hold on
plot(LT(v:w)/3600,A10(2,:),'b','LineWidth',2)   % Y rot leg1
plot(LT(v:w)/3600,A10(3,:),'g','LineWidth',2)   % X rot leg1
xlabel('Time, hrs')
ylabel('Angle, Degrees')
% title('Variation of the spherical joint angles w.r.t. time for leg 1')
legend({'Rot abt Y','Rot abt X'},'interpreter','latex')
grid on
% 

figure(16)
% subplot(2,3,2); 
% plot(LT(v:w)/3600,A20(1,:),'r*','LineWidth',2)       % Z rot leg2
hold on
plot(LT(v:w)/3600,A20(2,:),'b','LineWidth',2)       % Y rot leg2
plot(LT(v:w)/3600,A20(3,:),'g','LineWidth',2)       % X rot leg2
xlabel('Time, hrs')
ylabel('Angle, Degrees')
% title('Variation of the spherical joint angles w.r.t. time for leg 2')
legend({'Rot abt Y','Rot abt X'},'interpreter','latex')
grid on
% 

figure(17)
hold on
plot(LT(v:w)/3600,A30(2,:),'b','LineWidth',2)       % Y rot leg3
plot(LT(v:w)/3600,A30(3,:),'g','LineWidth',2)       % X rot leg3
xlabel('Time, hrs')
ylabel('Angle, Degrees')
% title('Variation of the spherical joint angles w.r.t. time for leg 3')
legend({'Rot abt Y','Rot abt X'},'interpreter','latex')
grid on
 
