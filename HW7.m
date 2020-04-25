%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HW7 Quick return mechanism
% Due date: 28-4-2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

TMT_method() % Function which computes all the values

w = 75;

theta2 = 0;
theta2d = w/60*2*pi;

start_angels;

theta4  = angels(1);
theta4d = angels(2);
theta5  = angels(3);
theta5d = angels(4);

yi = [theta2; theta2d; theta4; theta4d; theta5; theta5d];
n = 8;
t0 = 0;
tfinal = 1.8; % seconds
dt = tfinal/2^n;
yd = compute_yd(yi);
idx = 0;
time = 0:dt:tfinal;

%% RK4 method with coordinate projection
for idx = 1:length(time)
%calculate acceleration of angles, lambda's and slider
[k1,l1(:,idx),s1(:,idx)] = compute_yd(yi(:,idx));
[k2] = compute_yd(yi(:,idx) + dt/2*k1);
[k3] = compute_yd(yi(:,idx) + dt/2*k2);
[k4] = compute_yd(yi(:,idx) + dt*k3);
ydi(:,idx) = 1/6*(k1+2*k2+2*k3+k4);
yi(:,idx+1) = yi(:,idx) + dt*ydi(:,idx);

% Coordinate Projection Method

O2A  = 0.2; 
O4B  = 0.7; 
BC   = 0.6;
O4O2 = 0.3; 
O4G4 = 0.4; 
BG5  = 0.3; 
yC   = 0.9;

% Constrains from TMT_method copy
C1 = yC - O4B*sin(yi(3,idx))-BC*sin(yi(5,idx));
C2 = O2A*cos(yi(1,idx)) - sqrt(((O2A*cos(yi(1,idx)))^2 + (O4O2+O2A*sin(yi(1,idx)))^2))*cos(yi(3,idx));
Constraints = [C1;C2];

%jacobian copy of constrains
CCd = [0,-(7*cos(yi(3,idx)))/10, -(3*cos(yi(5,idx)))/5; (cos(yi(3,idx))*((2*cos(yi(1,idx))*sin(yi(1,idx)))/25 - (2*cos(yi(1,idx+1))*(sin(yi(1,idx))/5 + 2/5))/5))/(2*(cos(yi(1,idx))^2/25 + (sin(yi(1,idx))/5 + 2/5)^2)^(1/2)) - sin(yi(1,idx))/5, sin(yi(3,idx))*(cos(yi(1,idx))^2/25 + (sin(yi(1,idx))/5 + 2/5)^2)^(1/2), 0];
Ce = CCd.'*inv(CCd*CCd.'); % making C plus
error(:,idx) = Ce*-Constraints; % calculate the delta q

%calculate constrain forces formula 5.20
Cforce(:,idx) = -1*l1(:,idx).'*CCd;

% Adding delta q
yi(1,idx) = yi(1,idx)+error(1);
yi(3,idx) = yi(3,idx)+error(2);
yi(5,idx) = yi(5,idx)+error(3);
end




%% plotting system

% Question B, angluar speed of cranks
figure(1);
plot(time,ydi(1,:)); hold on
plot(time,ydi(3,:));
plot(time,ydi(5,:));
title(' Angle speed ')
xlabel('Time [sec] ')
ylabel('Speed rad/s')
legend('\theta 2', '\theta 3/4','\theta 5')

% Question C, sliding speed
figure(2);
subplot(2,2,1);
plot(time,s1(1,:)); hold on
title(' Slider 3 & rocker 4 ')
xlabel('Time [sec] ')
ylabel('Speed rad/s')
subplot(2,2,2);
plot(time,s1(2,:));
title(' Horizontal position slider 6 ')
xlabel('Time [sec] ')
ylabel('Position [m]')
subplot(2,2,3);
plot(time,s1(3,:));
title(' Speed slider 6')
xlabel('Time [sec] ')
ylabel('Speed [m/s]')
subplot(2,2,4);
plot(time,s1(4,:));
title(' Acceleration slider 6 ')
xlabel('Time [sec] ')
ylabel('Acceleration [m/s^2]')

% Question D, force plot
figure(3); 
plot(time,(Cforce(1,:))); hold on
plot(time,Cforce(3,:));
title(' Normal forces ')
xlabel('Time [sec] ')
ylabel('Force [N]')
legend('Slider 3 on 4','Slider 6')

% % simulation pot ctr+t to use
% Use after first run, running while
% 
% figure(4);
% axis equal
% hold on;
% pause(1)
% O2A = 0.2; O4B = 0.7; BC = 0.6;
% O4O2 = 0.3; O4G4 = 0.4; BG5 = 0.3; yC = 0.9;
% m3= 0.5; m4= 6; m5= 4; m6= 2; J4=10; J5= 6; F=1000;
% T=0; J2= 100; w= 75*2*pi/60; J3=0;
% p0 = plot(0,0,'o'); hold on;
% p2 = plot(0,O4O2,'o');



% for n = 1:length(yi)
% xA = O2A*cos(yi(1,n));
% yA = O4O2+O2A*sin(yi(1,n));
% xG4= O4G4*cos(yi(3,n));
% 
% yG4= O4G4*sin(yi(3,n));
% xB = O4B*cos(yi(3,n));
% yB = O4B*sin(yi(3,n));
% xG5= O4B*cos(yi(3,n))+BG5*cos(yi(5,n));
% yG5= O4B*sin(yi(3,n))+BG5*sin(yi(5,n));
% xC = O4B*cos(yi(3,n))+BC*cos(yi(5,n));
% yC = O4B*sin(yi(3,n))+BC*sin(yi(5,n));
% p3 = plot([0 xA],[O4O2 yA],'b-');
% p4 = plot([0 xB],[0 yB], 'b-');
% p5 = plot([xB xC],[yB yC],'b-'); 
% % pause time speed of simulation
% pause(0.001)
% delete(p3)
% delete(p4)
% delete(p5)
% xlim([-1 0.6])
% ylim([-0.5,1.5])
% end
function [yd, labda, Slider3] = compute_yd(yi)
% yi, current value -> to yd velocity and acc
theta2  = yi(1);
theta2d = yi(2);
theta4  = yi(3);
theta4d = yi(4);
theta5  = yi(5);
theta5d = yi(6);
  
symb_Acc ;
theta2dd = Acc(1);
theta4dd = Acc(2);
theta5dd = Acc(3);
labda = Acc(4:5,1);
symb_slid;
Slider3 = slid;
%projection version with constrain
yd = [theta2d; theta2dd; theta4d; theta4dd; theta5d; theta5dd];
end

%TMT method
function TMT_method()
syms theta2 theta4 theta5
syms theta2d theta4d theta5d theta2dd theta4dd theta5dd % for dx/dt I use xd etc
O2A  = 0.2; 
O4B  = 0.7; 
BC   = 0.6;
O4O2 = 0.3; 
O4G4 = 0.4; 
BG5  = 0.3; 
yC   = 0.9;

m3 = 0.5; 
m4 = 6;
m5 = 4; 
m6 = 2; 

J2 = 100;
J3 = 0;
J4 = 10; 
J5 = 6; 

F = 1000;
T = 0; 
 
w = 75*2*pi/60; 

%setting up kinematics.
q   = [theta2;theta4;theta5];
qd  = [theta2d;theta4d;theta5d];
qdd = [theta2dd; theta4dd; theta5dd];

xA = O2A*cos(theta2);
yA = O4O2+O2A*sin(theta2);
xG4= O4G4*cos(theta4);
yG4= O4G4*sin(theta4);
xG5= O4B*cos(theta4)+BG5*cos(theta5);
yG5= O4B*sin(theta4)+BG5*sin(theta5);
xC = O4B*cos(theta4)+BC*cos(theta5);

%making Ti,k matrix
Ki = [theta2;xA;yA;theta4;xG4;yG4;theta4;xG5;yG5;theta5;xC];%all kinematics combined
Ti = jacobian(Ki, q);
Tiv = Ti*qd; %velocity
gk = jacobian(Tiv,q)*qd;
Tacc = jacobian(Tiv,qd)*qdd + jacobian(Tiv,q)*qd;
%making Mass and gravity matrix and vector
M = diag([J2 m3 m3 J3 m4 m4 J4 m5 m5 J5 m6]);
% Constraints
C1 = yC - O4B*sin(theta4)-BC*sin(theta5);
C2 = O2A*cos(theta2)-sqrt((xA^2+yA^2))*cos(theta4);
Constraints = [C1;C2];
CCd = jacobian(Constraints,q);
CCdd = jacobian(CCd*qd,q)*qd;

%Forces
fi = [T; 0; 0; 0; 0; 0; 0; 0; 0; 0; F;];

%making TMT matrix
TMT = simplify([Ti.'*M*Ti CCd.';CCd zeros(2,2)]);
% dTdot = TMT
Q = Ti.'*(fi-M*gk);
Ftot = simplify([Q; -CCdd]);
Acc = TMT\Ftot;
if exist('symb_Acc.m', 'file')
! del symb_Acc.m
end
diary symb_Acc.m
disp('Acc = ['), disp(Acc), disp('];');
diary off
% making starting position with the first values.
theta4_start = atan((O4O2 + sin(theta2)*O2A)/(cos(theta2)*O2A));
theta4d_start = jacobian(theta4_start,q)*qd;
theta5_start = pi-asin((yC-O4B*sin(theta4_start))/BC);
theta5d_start = jacobian(theta5_start,q)*qd;
angels_start = [theta4_start; theta4d_start; theta5_start; theta5d_start];
if exist('start_angels.m', 'file')
! del start_angels.m
end
diary start_angels.m
disp('angels = ['), disp(angels_start), disp('];');
diary off
% Calculating speed of slider
slid(1) = sqrt(Tiv(2)^2+Tiv(3)^2); %sliding speed of slider 3
slid(2) = Ki(11);
slid(3) = Tiv(11);
slid(4) = Tacc(11);

if exist('symb_slid.m', 'file')
! del symb_slid.m
end
diary symb_slid.m
disp('slid = ['), disp(slid), disp('];');
diary off
end
