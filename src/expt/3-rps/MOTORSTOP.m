function MOTORSTOP
global ard
% global fileID
% fileID = fopen('EXP.txt','w');
% disp('ALL MOTORS STOPPED')
% fprintf(fileID,'ALL MOTORS STOPPED\n');

    % to stop leg 1
    ard.digitalWrite(30,0);
    ard.analogWrite(4,0);
    ard.digitalWrite(31,0);
    ard.analogWrite(5,0);

    % to stop leg 3
    ard.analogWrite(8,0);
    ard.analogWrite(9,0);
    ard.digitalWrite(34,0);
    ard.digitalWrite(35,0);
% 
    % to stop leg2 
    ard.analogWrite(6,0);
    ard.analogWrite(7,0);
    ard.digitalWrite(32,0);
    ard.digitalWrite(33,0);