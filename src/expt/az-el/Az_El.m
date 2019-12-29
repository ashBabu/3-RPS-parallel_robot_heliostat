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

year=2013;
month=3;
day=21;
d=date(year,month,day);

time_incre=900; % in seconds or in increments of 15 minutes
LT=0:time_incre:24*3600;

[elevation, azimuth]=alpha_azimuthcalc(LT,d,longitude,latitude);
[xs ,ys, zs]=sunposition(elevation, azimuth);  % xs,ys,zs are coordinates of the sun
OS=[xs;ys;zs];


% The lab wall is in the -ve x direction. So chosing a sun-vector in the
% evening, say, 3:15pm, to 4:15. Also, the portion of
% lab wall corresponding to the above time steps is neat and clean.
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
% 37-67 % 9am - 4.30pm

% To mark on the lab roof which is at 4.2m from the chosen global
% origin, the OS vector has to be multiplied by a constant so as to get the
% z as 4.2m
v_900 = OS(:,t_900)*(4.20/OS(3,t_900));
v_1000 = OS(:,t_1000)*(4.20/OS(3,t_1000));
v_1100 = OS(:,t_1100)*(4.20/OS(3,t_1100));
v_1200 = OS(:,t_1200)*(4.20/OS(3,t_1200));
v_1300 = OS(:,t_1300)*(4.20/OS(3,t_1300));

v_200 = OS(:,t_200)*(4.20/OS(3,t_200));
v_215 = OS(:,t_215)*(4.20/OS(3,t_215));
v_230 = OS(:,t_230)*(4.20/OS(3,t_230));
v_245 = OS(:,t_245)*(4.20/OS(3,t_245));
v_300 = OS(:,t_300)*(4.20/OS(3,t_300));

% To mark on the lab wall which is at -3.545m from the chosen global
% origin, the OS vector has to be multiplied by a constant so as to get the
% x as -3.545m
v_315 = OS(:,t_315)*(-3.540/OS(1,t_315));
v_330 = OS(:,t_330)*(-3.540/OS(1,t_330));
v_345 = OS(:,t_345)*(-3.540/OS(1,t_345));
v_400 = OS(:,t_400)*(-3.540/OS(1,t_400));
v_415 = OS(:,t_415)*(-3.540/OS(1,t_415));
v_430 = OS(:,t_430)*(-3.540 /OS(1,t_430));

% Tgt =[v_900 v_1000 v_1100 v_1200 v_1300 v_215 v_230 v_245 v_300 v_315 v_330 v_345 v_400 v_415]

%%
Motor_angle_Az_El =[]; Az_El_count =[];
% xs = 5.55;ys= -1.336; zs = 2.62;  % LASER coordinates
xl = 6.38;  yl = -1.21;  zl = 2.56;  % laser coordinates mounted on the wall
zG = 1.58;
B0t = [0;0;zG];   % Heliostat centre. 
zz = t_900;

for i=1:length(t_900:t_530)
    GP = OS(:,zz);
%     OP = OS(:,66);
%     GP=(OP-B0t)/norm((OP-B0t),2);
    
%     GP = OS(:,66);
    OL=[xl;yl;zl];  % laser vector
    GS=(OL-B0t)/norm((OL-B0t),2);
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
    Motor_angle_azimuth = N_azimuth  ;
    Motor_angle_elevation = 90-N_elevation ;
    

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

for i=1:length(t_900:t_530);
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
    pause(5)
    %    enc_A(i)= (str2double(g.command('TPA')))/9091 ;
    %    enc_B(i)= (str2double(g.command('TPB')))/9091 ;
end

g.command('PAA=0');
g.command('PAB=0');
g.command('BGA');
g.command('BGB');

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

figure(2)
plot(LT(t_900:t_415)/3600,Motor_angle_Az_El(:,2),'r','LineWidth',2)
hold on
plot(LT(t_900:t_415)/3600,Motor_angle_Az_El(:,1),'b','LineWidth',2)
plot(LT(t_900:t_415)/3600,-enc_B,'r*','LineWidth',2)
plot(LT(t_900:t_415)/3600,enc_A,'b*','LineWidth',2)
xlabel('Time  \it{t / hr}')
ylabel('Motor angles  \it{\theta / (^0)}')
title('Actuation required for the motor')
lgnd=legend('Ref. Elevation',' Ref. Azimuth','Encoder Elevation','Encoder Azimuth');
grid on

g.command('TPB=0')

g.command('SH')
g.command('PRB=80000')
g.command('BGB')

g.command('MOB')

enc_counts = [40000 80000 120000 160000 200000 240000 200000 160000 120000 80000 30000 0 30000 80000 30000 110000]';
acc = [0.3 9.8; 1.2 9.8; 1.9 9.6; 2.8 9.5; 3.6 9.2; 4.3 8.9; 4.9 8.6; 4.3 8.9; 3.6 9.2; 2.9 9.5; 2.1 9.6; 1.2 9.7; 0.7 9.8;1.2 9.8; 2.2 9.6; 1.3 9.8; 2.8 9.5];

% aas = [0,0;5.22765254675825,40000;9.44170656622844,80000;14.6687826382438,120000;19.6172174092716,160000;24.0339233227428,200000;27.9196772520063,240000;24.0339233227428,200000;19.6172174092716,160000;15.2220946078582,120000;10.5856824182547,80000;5.29890114776965,30000;0.2,0;5.22765254675825,30000;11.1540038111943,80000;5.80293919636792,30000;14.6687826382438,110000]
jj = 1;
theta = [];
for i=1: length(acc)
   ang(i) = atan2d(acc(i,1),acc(i,2)) ;
   if i>1
   theta = cat(1,theta,ang(i)-ang(1));
   
   end
end

 plot(theta,enc_counts, 'r*')


clc
 ff=fit(theta,enc_counts,'poly1')
 ff.p1
 ff.p2
