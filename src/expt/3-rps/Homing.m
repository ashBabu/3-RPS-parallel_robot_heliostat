% delete(instrfind({'PORT'},{'COM7'}))
% a = arduino('COM7');
Act_tme = [0 -1 1];
tic
while(toc <= 20)
    %% Leg1
    if(Act_tme(1) <0)
        % bck
        ard.digitalWrite(30,0);
        ard.analogWrite(4,0);
        ard.digitalWrite(31,1);
        ard.analogWrite(5,200);
    %
    else
        % fwd
        ard.digitalWrite(31,0);
        ard.analogWrite(5,0);
        ard.digitalWrite(30,1);
        ard.analogWrite(4,200);
    end
    %% Leg 3
    if(Act_tme(3) <0)
    % bwd
        ard.digitalWrite(34,0);
        ard.analogWrite(8,0);
        ard.digitalWrite(35,1);
        ard.analogWrite(9,200);
    else
%     % fwd
    ard.digitalWrite(35,0);
    ard.analogWrite(9,0);
    ard.digitalWrite(34,1);
    ard.analogWrite(8,200);
    
    end
     if(Act_tme(2) <0)
       %% Leg 2
    
    % bwd
        ard.digitalWrite(32,0);
        ard.analogWrite(6,0);
        ard.digitalWrite(33,1);
        ard.analogWrite(7,200);
    else
%     % fwd
    ard.digitalWrite(33,0);
    ard.analogWrite(7,0);
    ard.digitalWrite(32,1);
    ard.analogWrite(6,200);
    
    end
end
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

% to stop leg2 
ard.analogWrite(6,0);
ard.analogWrite(7,0);
ard.digitalWrite(32,0);
ard.digitalWrite(33,0);