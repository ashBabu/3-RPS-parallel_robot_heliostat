%inputs needed are station co-ordinates longitude,latitude in degrees,min,sec ...
%altitude in kms and time of proposed sun visibility in year,month,day in Christian era,
%hr,min and sec in GMT (NOT BASED ON LOCAL TIME).
%Input the following
clear all
home
close all
clc
%% The following part of the code is used to mark point on the lab wall to show sun vector
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
month=11;
day=10;
d=date(year,month,day);

time_incre=900; % in seconds or in increments of 15 minutes
LT=0:time_incre:24*3600;

[elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude);
[xs ,ys, zs]=sunposition(elevation, azimuth);  % xs,ys,zs are coordinates of the sun
OS=[xs;ys;zs];

t_900 = 9/time_incre*3600+1;
t_915 = t_900 + 1;
t_930 = t_900 + 2;
t_945 = t_900 + 3;
t_1000 = 10/time_incre*3600+1;
t_1015 = t_1000 + 1;
t_1030 = t_1000 + 2;
t_1045 = t_1000 + 3;

t_1100 = 11/time_incre*3600+1;
t_1115 = t_1100 + 1;
t_1130 = t_1100 + 2;
t_1145 = t_1100 + 3;

t_1200 = 12/time_incre*3600+1;
t_1215 = t_1200 + 1;
t_1230 = t_1200 + 2;
t_1245 = t_1200 + 3;

t_1300 = 13/time_incre*3600+1;
t_1315 = t_1300 + 1;
t_1330 = t_1300 + 2;
t_1345 = t_1300 + 3;

t_200 = 14/time_incre*3600+1;
t_215 = t_200 + 1;
t_230 = t_200 + 2;
t_245 = t_200 + 3;
t_300 = t_200 + 4;
t_315 = t_200 + 5;
t_330 = t_200 + 6;
t_345 = t_200 + 7;
t_400 = 16/time_incre*3600+1;
t_415 = t_400 + 1;
t_430 = t_400 + 2;
t_445 = t_400 + 3;
t_500 = 17/time_incre*3600+1;
t_515 = t_500 + 1;
t_530 = t_500 + 2;
%% Plot of azimuth and elevation angles
figure(1)
plot(LT(t_900:t_415)/3600,elevation((t_900:t_415)),'r','LineWidth',4)
hold on
plot(LT(t_900:t_415)/3600,azimuth((t_900:t_415)),'b','LineWidth',4)
xlabel('Time  \it{t / hr}')
ylabel('Sun vector angles  \it{\theta / (^0)}')
title('Azimuth and Elevation angles of the sun')
lgnd=legend('Elevation',' Azimuth');
grid on
%%
Motor_angle_Az_El =[]; Az_El_count =[];
zG = 1.58;
zz = t_900;

a = 0;  b = 0;  c = 6.72; % wrt TCS
OP = [a;b;c];  % wrt TCS; O is the origin of TCS which is R[10deg, Z] wrt gcs.
OO1 = [-14;5.45;0];
O1G = [0;0;zG];   % Heliostat centre. same in both TCS,gcs,bcs
OG = (OO1+O1G);
GP =(OP-OG)/norm((OP-OG),2);
% Rot = eye(3);
incl = 10;  % angle by which TCS and gcs are rotated
Rot = [cosd(incl) -sind(incl) 0;
       sind(incl)  cosd(incl) 0;
           0           0      1];
        
for i=1:length(t_1200:t_530)
    GS=Rot'*OS(:,zz);
    n = GS + GP;
    normal=n/norm(n,2);
    xn=normal(1);   yn=normal(2);       zn=normal(3); qq=sqrt(xn^2+yn^2);
    rotmat =[xn*zn/qq -yn/qq xn ;yn*zn/qq xn/qq  yn;-(xn^2+yn^2)/qq 0 zn];
        
    %% To find the azimuth and elevation angles of the heliostat normal
    
    N_azimuth= atan2d(normal(2),normal(1));
    l=sqrt(normal(2)^2+normal(1)^2);
    N_elevation = atand(normal(3)/l);  % ANGLE OF THE MIRRO NORMAL WRT HORIZONTAL.
                            %IN THE AZ-EL HELIOSTAT MADE IN THE LAB, THE MOTOR ANGLE
                            % SHOULD BE 90-N_elevation
    
%     Motor_angle_azimuth = N_azimuth - 3 ;
%     Motor_angle_elevation = 90-N_elevation - 6;
%     Motor_angle_azimuth = N_azimuth - 3 ;
%     Motor_angle_elevation = 90-N_elevation - 5;
    Motor_angle_azimuth = N_azimuth-10 ;
    Motor_angle_elevation = 90-N_elevation-18 ;
    

    counts_per_degree = 9091;
    % 1 degree = 9091 counts of the motor
    azi_count  = counts_per_degree*Motor_angle_azimuth;
    ele_count  = counts_per_degree*Motor_angle_elevation;
    
    Motor_angle_Az_El  = cat(1,Motor_angle_Az_El,[Motor_angle_azimuth Motor_angle_elevation])
    Az_El_count = cat(1, Az_El_count, [azi_count -ele_count])  % elevation is multiplied by a -ve sign
                                            % bcoz the Az-El heliostat made
                                            % in the lab rotates in -ve
                                            % direction if a positive value
                                            % is given and vice versa. This
                                            % is bcoz of the direction of
                                            % the axis chosen.
    zz=zz+1;
end


% Az_El_count = [100000 100000; -100000 -100000]
%% for Galil Motion Control
g = actxserver('galil');%set the variable g to the GalilTools COM wrapper 
response = g.libraryVersion;%Retrieve the GalilTools library versions 
disp(response);%display GalilTools library version 
% g.address = '';%Open connections dialog box 
g.address = '192.168.1.1';%Open connections dialog box 

g.command('TPA=0');
g.command('TPB=0');
g.command('SH');


[srt_azel_cnt,Ind] = sort(abs(Az_El_count),2);

for i=1:length(t_1200:t_530);
    Enc_cnt = 0.00;
    str1 = ['''PAA=',num2str(Az_El_count(i,1)),''''];
    %     str2 = ['''PA',num2str(Az_El_count(i,:)),''''];
    str2 = ['''PAB=',num2str(Az_El_count(i,2)),''''];
    g.command(eval(str1))
    g.command(eval(str2))
    g.command('BGA');
    %     g.command('AMA');
    g.command('BGB');
    while (abs(Enc_cnt) <=abs(srt_azel_cnt(i,end)))
        if Ind(i,end) ==1
            Enc_cnt = str2num(g.command('TPA'))
        elseif Ind(i,end) == 2
            Enc_cnt = str2num(g.command('TPB'))
        end
        %             break;
        %             pause(1)
    end
    pause
    %    enc_A(i)= (str2double(g.command('TPA')))/9091 ;
    %    enc_B(i)= (str2double(g.command('TPB')))/9091 ;
end

% g.command('PAA=0');
% g.command('PAB=0');
% g.command('BGA');
% g.command('BGB');

% for i=1:length(Az_El_count)
    
%     str1 = ['''PRA=',num2str(Az_El_count(i,1)),''''];
% %     str2 = ['''PA',num2str(Az_El_count(i,:)),''''];
%     str2 = ['''PRB=',num2str(Az_El_count(i,2)),''''];
%     g.command(eval(str1))
%     g.command(eval(str2))
%     g.command('BGA');
%     while (str2double(Enc_cnt_A) <=Az_El_count(i,1))
%         
% %         g.command('AMA');
%         Enc_cnt_A = g.command('TPA')
%         break;
% %         pause(1)
%     end
%     tt=tt+1
% %     while (str2double(g.command('TPB=?')) <=Az_El_count(i,2))
% %         g.command('BGB');
% % %         g.command('AMB');
% %         pause(1)
% %     end
%     pause(6)
%     
% end
% g.command('ST');
% g.command('MO');
% i=1;


% % g.command('MO');

% figure(2)
% plot(LT(t_900:t_415)/3600,Motor_angle_Az_El(:,2),'r','LineWidth',2)
% hold on
% plot(LT(t_900:t_415)/3600,Motor_angle_Az_El(:,1),'b','LineWidth',2)
% plot(LT(t_900:t_415)/3600,-enc_B,'r*','LineWidth',2)
% plot(LT(t_900:t_415)/3600,enc_A,'b*','LineWidth',2)
% xlabel('Time  \it{t / hr}')
% ylabel('Motor angles  \it{\theta / (^0)}')
% title('Actuation required for the motor')
% lgnd=legend('Ref. Elevation',' Ref. Azimuth','Encoder Elevation','Encoder Azimuth');
% grid on

% g.command('TPB=0')
% 
% g.command('SH')
% g.command('PRB=80000')
% g.command('BGB')
% 
% g.command('MOB')
% 
% enc_counts = [40000 80000 120000 160000 200000 240000 200000 160000 120000 80000 30000 0 30000 80000 30000 110000]';
% acc = [0.3 9.8; 1.2 9.8; 1.9 9.6; 2.8 9.5; 3.6 9.2; 4.3 8.9; 4.9 8.6; 4.3 8.9; 3.6 9.2; 2.9 9.5; 2.1 9.6; 1.2 9.7; 0.7 9.8;1.2 9.8; 2.2 9.6; 1.3 9.8; 2.8 9.5];
% 
% % aas = [0,0;5.22765254675825,40000;9.44170656622844,80000;14.6687826382438,120000;19.6172174092716,160000;24.0339233227428,200000;27.9196772520063,240000;24.0339233227428,200000;19.6172174092716,160000;15.2220946078582,120000;10.5856824182547,80000;5.29890114776965,30000;0.2,0;5.22765254675825,30000;11.1540038111943,80000;5.80293919636792,30000;14.6687826382438,110000]
% jj = 1;
% theta = [];
% for i=1: length(acc)
%    ang(i) = atan2d(acc(i,1),acc(i,2)) ;
%    if i>1
%    theta = cat(1,theta,ang(i)-ang(1));
%    
%    end
% end
% 
%  plot(theta,enc_counts, 'r*')
% 
% 
% clc
%  ff=fit(theta,enc_counts,'poly1')
%  ff.p1
%  ff.p2
