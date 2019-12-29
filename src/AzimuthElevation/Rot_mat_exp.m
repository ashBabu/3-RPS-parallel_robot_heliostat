axis = [0;0;1];
Skew_axis_matrix = [0 -axis(3) axis(2);
                    axis(3) 0 -axis(1);
                    -axis(2) axis(1) 0];
 Rz = expm(180*Skew_axis_matrix)               