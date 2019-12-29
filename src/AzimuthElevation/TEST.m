% This program uses fsolve to solve for AZ-EL case. Since there are only
% two rotations happening, viz. Rz and Ry, in AZ-EL case, the R(3,2) element of the rotation
% matrix remains zero. The sun vector, the reflected ray are completely
% known(since the center is fixed). Hence, the direction cosines of the
% normal can be found out. Therefore there are only 5 unknowns, viz.,
% n1,n2,n3,o1 and o2 which are related by 5 equations, viz., n.n=1, o.o=1,
% n.a=0, n.o=0, a.o=0 ; Fsolve is used to solve this. It should be noted
% that depending upon the values of the initial condition, x0, the solution
% changes 
clear all
home
close all
global O1M OP Rot_gam OO1
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

a=0;b=0;c=65;  % Receiver coordinates wrt gcs
OP = [a;b;c];  % Receiver vector
year=2016;
month=11;
day=1;
d=date(year,month,day);

time_incre=10*60;  % seconds
LT=0*3600:time_incre:24*3600;

[elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude);  % aximuth and elevation of sun
[xs, ys, zs]=sunposition(elevation, azimuth);             % x2,y2,z2 are coordinates of the sun
OS=[xs;ys;zs];  % Sun vector
rad=100; % radius of the mirror from the origin
psi=30; % psi is the angle the heliostat centre makes with the local east(X) axis
zG=.5;  % height of the heliostat

KJ = -.25; % distance between the mirror center and the top platform center

O1M = [0;0;zG];   % vector from the base origin to the mirror origin described wrt bcs
OO1=[rad*cosd(psi);rad*sind(psi);0];   % vector from the global origin to the base origin
v=9/time_incre*3600+1; % start time for the plots
w=16/time_incre*3600+1; % end time for the plots
o3=0;  % three-two element of the rotation matrix

gamma = 0; % the angle by which the x_b is inclined towards X (East)
Rot_gam=[ cosd(gamma) -sind(gamma) 0;
          sind(gamma)  cosd(gamma) 0;
              0            0       1]; 

% x0 = [0.2;0.61;0.5;-0.01;0.015]; 
x0 = [.7;-.6;0.1;-0.7;-0.01;-0.05];
% x0 = [.7;-.6;0.1;-0.7;-0.01;-0.05];
% x0 = [0.2903; -0.4146; -0.8625; 0.8192; 0.5735];  % for 3rd quadrant
% x0 = [0.2903; -0.4146; 0.5625; -0.8192; 0.5735];  
% x0 = [0.1761;-0.5876;-0.7897;0.9579;0.2871]; % Initial guesses for n1,n2,n3,o1 and o2
options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt');
for i=1:length(v:w)
    zz=v-1+i;
    fhandle=@(x)for_mainprog_AZ_EL_3_UPU(x,OS(:,zz));
    [x,fval] = fsolve(fhandle,x0,options);
    x0=x;
    n1=x(1);n2=x(2);n3=x(3);o1=x(4);o2=x(5);
    O1P_bcs=Rot_gam'*(OP-OO1);      % O1P_bcs = O1P described in bcs
    g=O1P_bcs-O1M;   % wrt bcs
    GP=g/norm(g,2);  % unit reflected ray wrt bcs
    GS=Rot_gam'*OS(:,zz);         % Sun vector in bcs
    GN=(GS+GP)/norm((GS+GP),2);   % unit mirror normal wrt bcs  
    a1=GN(1); a2=GN(2); a3=GN(3); 

    Trans_matr(:,:,i)= [n1 o1 a1 O1M(1);      % transformation which takes mirror to base coordinate system
                        n2 o2 a2 O1M(2);
                        n3 o3 a3 O1M(3);
                         0   0  0   1] ;
    Rot(:,:,i)= [n1 o1 a1 ;              % rotation which takes mirror to base coordinate system
                 n2 o2 a2 ;
                 n3 o3 a3];  
%             R(:,:,i) = Rot_gam*Rot(:,:,i);   % rotation which takes mirror to global coordinate system
            O1G(:,i) = O1M+KJ*Rot(:,3,i);   % Vector from the base origin to the origin of the top platform
%             plot3(O1G(1,i),O1G(2,i),O1G(3,i),'r*')
%             hold on
end

%% Leg lengths
figure(1)
rb=.5;  % circum radius of the base equilateral triangle, B1B2/sqrt(3) in meters
rp=0.25;  % circum radius of the platform,P1P2/sqrt(3)in meters

ss=0.08; % square mirror half length
% M is the centre of the mirror and M1, M2, M3, and M4 are the corners
MM1=[ss;ss;0];  MM2=[-ss;ss;0];MM3=[-ss;-ss;0];  MM4=[ss;-ss;0];

O1R1=[rb;0;0];O1R2=[-0.5*rb;sqrt(3)*rb/2;0];O1R3=[-0.5*rb;-sqrt(3)*rb/2;0];
O1Rs=[O1R1 O1R2 O1R3];

GS1=[rp;0;0];GS2=[-0.5*rp;sqrt(3)*rp/2;0];GS3=[-0.5*rp;-sqrt(3)*rp/2;0];
GSs=[GS1 GS2 GS3];
for i=1:length(Rot)
    
    O1S1=Rot(:,:,i)*GS1+O1G(:,i);   % w.r.t base coordinate system
    l1_dir(:,i)=O1S1-O1R1;
    l1(i,:)=norm(l1_dir(:,i));
        
    O1S2=Rot(:,:,i)*GS2+O1G(:,i);
    l2_dir(:,i)=O1S2-O1R2;
    l2(i,:)=norm(l2_dir(:,i));
        
    O1S3=Rot(:,:,i)*GS3+O1G(:,i);
    l3_dir(:,i)=O1S3-O1R3;
    l3(i,:)=norm(l3_dir(:,i));
end

plot(LT(v:w)/3600, l1*1000,'r','LineWidth', 2)
hold on
plot(LT(v:w)/3600, l2*1000,'g','LineWidth', 2)
plot(LT(v:w)/3600, l3*1000,'b','LineWidth', 2)
xlabel('Time, hrs')
ylabel('leg length, mm')
% title('leg length Vs Time Az-El of 3-UPU')
legend('leg1','leg2', 'leg3')

leg_vec = [LT(v:w)'/3600 l1*1000 l2*1000 l3*1000];
leg_vec(1:40,:) = []; 
load('Rot_mat_TA.mat')
for i=1:length(Rot)
   Rzz(:,:,i) = Rot(:,:,i)'*Rot_TA(:,:,i);
   angle(i,:) = atan2d(Rzz(2,1,i),Rzz(1,1,i));
end

%% For Animation
TransmatrN=[Rot_gam,OO1;[0 0 0 1]];           % transformation which takes base to global coordinate system
base=TransmatrN*[O1Rs;[1 1 1]];
writerObj=VideoWriter('3UPU.avi');
% writerObj.FrameRate = 9;
open(writerObj);
figure(2)
for i=1:length(Trans_matr)% 
% for i =1:5
    xxx=base(1,:);                                   % wrt global coordinate system
    yyy=base(2,:);                                   % wrt global coordinate system
    zzz=base(3,:);                                   % wrt global coordinate system
    fill3(xxx,yyy,zzz,'r')
    hold on
    OM = OO1 + Rot_gam*O1M;
    OM1 = OO1 + Rot_gam*(O1M+Rot(:,:,i)*MM1);
    OM2 = OO1 + Rot_gam*(O1M+Rot(:,:,i)*MM2);
    OM3 = OO1 + Rot_gam*(O1M+Rot(:,:,i)*MM3);
    OM4 = OO1 + Rot_gam*(O1M+Rot(:,:,i)*MM4);
    
    Mirror=[OM1 OM2 OM3 OM4];
    xm=Mirror(1,:);
    ym=Mirror(2,:);
    zm=Mirror(3,:);
%     fill3(xm,ym,zm,'g')  % top platform or mirror

    OS1=OO1+Rot_gam*(O1G(:,i)+Rot(:,:,i)*GS1);
    OS2=OO1+Rot_gam*(O1G(:,i)+Rot(:,:,i)*GS2);
    OS3=OO1+Rot_gam*(O1G(:,i)+Rot(:,:,i)*GS3);
    top_triangle=[OS1 OS2 OS3];
    xt=top_triangle(1,:);
    yt=top_triangle(2,:);
    zt=top_triangle(3,:);
    fill3(xt,yt,zt,'c')  % top triangle
    
    TRANSMATR=Rot_gam*Rot(:,:,i);
    plot3([xt(1) xxx(1)],[yt(1) yyy(1)],[zt(1) zzz(1)],'--k','LineWidth',2.5)
    plot3([xt(2) xxx(2)],[yt(2) yyy(2)],[zt(2) zzz(2)],'--k','LineWidth',2.5)
    plot3([xt(3) xxx(3)],[yt(3) yyy(3)],[zt(3) zzz(3)],'--k','LineWidth',2.5)
    %    %% changing wrt global coordinate system
    
    kons=1;kon=.5;zz=v-1+i;
    normal=plot3([OM(1,1) OM(1,1)+kon*TRANSMATR(1,3)],[OM(2,1) OM(2,1)+kon*TRANSMATR(2,3)],[OM(3,1) OM(3,1)+kon*TRANSMATR(3,3)],'r','LineWidth',2);
    %     proj=plot3([B0t(1,i) B0t(1,i)+kon*Trans_matr(1,3,i)],[B0t(2,i) B0t(2,i)+kon*Trans_matr(2,3,i)],[0 0],'r','LineWidth',2);
    
    sunposition=plot3([OM(1,1) OM(1,1)+kons*xs(zz)],[OM(2,1) OM(2,1)+kons*ys(zz)],[OM(3,1) OM(3,1)+kons*zs(zz)],'k','LineWidth',2);   % position of the sun at every time_incre looking from plate coordinte system
    plot3(OM(1,1)+kons*xs(v:zz),OM(2,1)+kons*ys(v:zz),OM(3,1)+kons*zs(v:zz),'c*');  % position of the sun at every time_incre looking from plate coordinte system
    
    %%
    kjs =.5;
    plot3([OM(1,1) OM(1,1)+kjs*GP(1,1)],[OM(2,1) OM(2,1)+kjs*GP(2,1)],[OM(3,1) OM(3,1)+kjs*GP(3,1)],'b','LineWidth',2);

    %%
    % %     coordintes of the receiver
%     plot3([OM(1,1) a],[OM(2,1) b],[OM(3,1) c],'b','LineWidth',2);
%     plot3([OM(1,1) a],[OM(2,1) b],[OM(3,1) c],'k*');
%     plot3([a a],[b b],[0 c],'b');
    
    %     xlim([10 30]); ylim([40 50]); zlim([0 20]);
    %     xlim([-50 120]); ylim([0 30]); zlim([0 50]);
    mg = .3;
        axis([OO1(1)-mg OO1(1)+mg OO1(2)-mg OO1(2)+mg 0 1]);  
%         axis([7 10 4.5 5.5 0 1.5]); 
%     view(82,6)            % to see the
    %     Receiver
    axis square
%     axis tight
%     view(-26,34);
    grid on
    xlabel('East(x)')
    ylabel('North(y)')
    zlabel('Zenith(p)')
    hold off
    M(i)=getframe(gcf);
        writeVideo(writerObj,M(i));
    %     pause(0.1)
end 
close(writerObj)