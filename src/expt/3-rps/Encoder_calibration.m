close all
clear all
% Encoder calibration
time_fwd = [0 3 6 9 13 15 17 19 22 27]';
time_bwd = [0 5 8 10 12 14 18 21 24 27]';

%%
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 15)
set(0,'DefaultAxesFontWeight', 'Bold')

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 13)

%% leg 1

Dist_fwd_l1 = [0 8 15 22 30 34 39 44 51 62]'; %  forward or outward motion of leg 1
Enc_cnt_fwd_l1 = [0 674 1344 2039 2937 3423 3877 4309 4949 6130]'; % Encoder counts of leg 1 in forward motion

Dist_bwd_l1 = [0 -12 -19 -24 -29 -33 -42 -49 -57 -64]'; %  backward or inward motion of leg 1
Enc_cnt_bwd_l1 = [0 -3006 -4779 -5921 -6996 -8143 -10507 -11677 -13364 -15032]'; % Encoder counts of leg 1 in bacward motion

l1_fwd = fit(time_fwd,Dist_fwd_l1, 'poly1');
l1_bwd = fit(time_bwd,Dist_bwd_l1, 'poly1');

figure(1)
subplot(1,2,1)
plot(l1_fwd,time_fwd,Dist_fwd_l1,'r*')
a = l1_fwd.p1;
b = l1_fwd.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(5,50,polyfit_str);
ylabel('Distance (mm)')
xlabel('time (s)')
title ('Forward motion for leg 1')

subplot(1,2,2)
plot(l1_bwd,time_bwd,Dist_bwd_l1,'b*')
a = l1_bwd.p1;
b = l1_bwd.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(5,-50,polyfit_str);
xlabel('time (s)')
ylabel('Distance (mm)')
title ('Backward motion  for leg 1')

% Dist vs Encoder counts leg 1

l1_fwd_en = fit(Dist_fwd_l1,Enc_cnt_fwd_l1, 'poly1');
l1_bwd_en = fit(Dist_bwd_l1,Enc_cnt_bwd_l1, 'poly1');

figure(2)
subplot(1,2,1)
plot(l1_fwd_en,Dist_fwd_l1,Enc_cnt_fwd_l1,'r*')
a = l1_fwd_en.p1;
b = l1_fwd_en.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(5,3000,polyfit_str);
xlabel('Distance (mm)')
ylabel('Encoder counts')
title ('Forward motion for leg 1')

subplot(1,2,2)
plot(l1_bwd_en,Dist_bwd_l1,Enc_cnt_bwd_l1,'b*')
a = l1_bwd_en.p1;
b = l1_bwd_en.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(-65,-6000,polyfit_str);
ylabel('Encoder counts')
xlabel('Distance (mm)')
title ('Backward motion  for leg 1')


%% leg 2

Dist_fwd_l2 = [0 8 15 22 31 35 40 45 51 62]'; %  forward or outward motion of leg 1
Enc_cnt_fwd_l2 = [0 5225 10434 15584 22468 25810 29187 32489 37697 46377]'; % Encoder counts of leg 1 in forward motion

Dist_bwd_l2 = [0 -11 -19 -23 -29 -34 -43 -49 -56 -63]'; %  backward or inward motion of leg 1
Enc_cnt_bwd_l2 = [0 -8989 -14366 -17885 -21349 -24934 -32256 -37652 -43102 -48522]'; % Encoder counts of leg 1 in bacward motion

l2_fwd = fit(time_fwd,Dist_fwd_l2, 'poly1');
l2_bwd = fit(time_bwd,Dist_bwd_l2, 'poly1');

figure(3)
subplot(1,2,1)
plot(l2_fwd,time_fwd,Dist_fwd_l2,'r*')
a = l2_fwd.p1;
b = l2_fwd.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(5,50,polyfit_str);
ylabel('Distance (mm)')
xlabel('time (s)')
title ('Forward motion for leg 2')

subplot(1,2,2)
plot(l2_bwd,time_bwd,Dist_bwd_l2,'b*')
a = l2_bwd.p1;
b = l2_bwd.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(5,-50,polyfit_str);
xlabel('time (s)')
ylabel('Distance (mm)')
title ('Backward motion  for leg 2')

% Dist vs Encoder counts leg 2

l2_fwd_en = fit(Dist_fwd_l2,Enc_cnt_fwd_l2, 'poly1');
l2_bwd_en = fit(Dist_bwd_l2,Enc_cnt_bwd_l2, 'poly1');

figure(4)
subplot(1,2,1)
plot(l2_fwd_en,Dist_fwd_l2,Enc_cnt_fwd_l2,'r*')
a = l2_fwd_en.p1;
b = l2_fwd_en.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(5,3000,polyfit_str);
xlabel('Distance (mm)')
ylabel('Encoder counts')
title ('Forward motion for leg 2')

subplot(1,2,2)
plot(l2_bwd_en,Dist_bwd_l2,Enc_cnt_bwd_l2,'b*')
a = l2_bwd_en.p1;
b = l2_bwd_en.p2;
polyfit_str = ['y = ' num2str(a) ' *x + ' num2str(b)];
% polyfit_str qill be : y = 4*x + 2
text(-65,-6000,polyfit_str);
ylabel('Encoder counts')
xlabel('Distance (mm)')
title ('Backward motion  for leg 2')



