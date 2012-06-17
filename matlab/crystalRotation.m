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

C2G*C2G'-eye(3)