clc

delta = [1.0 1.0 1.0]';

delta_AC_CB_BA = [0.0,-1.0, 1.0,	
                  1.0, 0.0,-1.0,	
                 -1.0, 1.0, 0.0];	
%Eigen::Matrix<double,dim,dim> R;
R(:,3) = delta/norm(delta);
R(:,1) = delta_AC_CB_BA(:,1)/norm(delta_AC_CB_BA(:,1));
R(:,2) = cross(R(:,3),R(:,1));

format long
C2G=R'


R(:,2) = delta/norm(delta);
R(:,1) = delta_AC_CB_BA(:,1)/norm(delta_AC_CB_BA(:,1));
R(:,3) = cross(R(:,1),R(:,2));

format long
C2G=R'
