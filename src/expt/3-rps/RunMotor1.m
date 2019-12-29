function RunMotor1(act_time)
global ard
% fileID = fopen('EXP.txt','w');
% global fileID
    if act_time <0
        % bck
%         disp('Motor 1 back')
%         fprintf(fileID,'Motor 1 back\n');
        ard.digitalWrite(30,0);
        ard.analogWrite(4,0);
        ard.digitalWrite(31,1);
        ard.analogWrite(5,200);
    else
        % fwd
%         disp('Motor 1 fwd')
%         fprintf(fileID,'Motor 1 fwd\n');
        ard.digitalWrite(31,0);
        ard.analogWrite(5,0);
        ard.digitalWrite(30,1);
        ard.analogWrite(4,200);
    end
    
