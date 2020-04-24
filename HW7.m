%% Assignment 7
clear all
close all
%functions are below.
%TMT calculates all values and write them to another file
%( symb_acc gives
%acc of angles, Symb_angels.m gives the starting position of
%all angles, symb_slid gives slider acceleration)
TMT_1()
g = 0;
w = 75;
gd = w/60*2*pi;
start_angels;
a = angels(1);
ad = angels(2);
b = angels(3);
bd = angels(4);
yi = [g; gd; a; ad; b; bd];
n =8;
t0 = 0;
tfinal = 1.8; % two rounds is about 1.7... sec so i take 1.8s
dt = tfinal/2^n;
ydot = calcydot(yi);
idx = 0;
time = 0:dt:tfinal;
%% RK4 method with coordinate projection
for idx = 1:length(time)
%calculate acceleration of angles, lambda's and slider
[k1,l1(:,idx),s1(:,idx)] = calcydot(yi(:,idx));
[k2] = calcydot(yi(:,idx) + dt/2*k1);
[k3] = calcydot(yi(:,idx) + dt/2*k2);
[k4] = calcydot(yi(:,idx) + dt*k3);
ydi(:,idx) = 1/6*(k1+2*k2+2*k3+k4);
yi(:,idx+1) = yi(:,idx) + dt*ydi(:,idx);
% yi(:,idx+1) = yi(:,idx) + dt*k1; Euler method
% Coordinate Projection Method

O2A = 0.2; O4B = 0.7; BC = 0.6;
O4O2 = 0.3; O4G4 = 0.4; BG5 = 0.3; yC = 0.9;
% Constrains from TMT_1 copy
Cy = yC - O4B*sin(yi(3,idx))-BC*sin(yi(5,idx));
ConA = O2A*cos(yi(1,idx)) - sqrt(((O2A*cos(yi(1,idx)))^2 + (O4O2+O2A*sin(yi(1,idx)))^2))*cos(yi(3,idx));
Constra = [Cy;ConA];
%jacobian copy of constrains
CCd = [0,-(7*cos(yi(3,idx)))/10, -(3*cos(yi(5,idx)))/5; (cos(yi(3,idx))*((2*cos(yi(1,idx))*sin(yi(1,idx)))/25 - (2*cos(yi(1,idx+1))*(sin(yi(1,idx))/5 + 2/5))/5))/(2*(cos(yi(1,idx))^2/25 + (sin(yi(1,idx))/5 + 2/5)^2)^(1/2)) - sin(yi(1,idx))/5, sin(yi(3,idx))*(cos(yi(1,idx))^2/25 + (sin(yi(1,idx))/5 + 2/5)^2)^(1/2), 0];
Ce = CCd.'*inv(CCd*CCd.'); % making C plus
error(:,idx) = Ce*-Constra; % calculate the delta q
%calculate constrain forces formula 5.20
Cforce(:,idx) = -1*l1(:,idx).'*CCd;
% Adding delta q
yi(1,idx) = yi(1,idx)+error(1);
yi(3,idx) = yi(3,idx)+error(2);
yi(5,idx) = yi(5,idx)+error(3);
end
%% plotting system
figure(1);% angluar speed of cranks
plot(time,ydi(1,:)); hold on
plot(time,ydi(3,:));
plot(time,ydi(5,:));
title(' Angle speed ')
xlabel('Time [sec] ')
ylabel('Speed rad/s')
legend('\theta 2', '\theta 3/4','\theta 5')
figure(2);% sliding speed question c
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
figure(3); % force plot question d
% plot(time,sqrt(Cforce(1,:).^2+Cforce(2,:).^2)); hold on
plot(time,(Cforce(1,:))); hold on
plot(time,Cforce(3,:));
% plot(time,Cforce(2,:));
title(' Normal forces ')
xlabel('Time [sec] ')
ylabel('Force [N]')
legend('Slider 3 on 4','Slider 6')
% plot(time,Cforce(1,:));
%% simulation pot ctr+t to use
% Use after first run, running while
%
figure(4);
axis equal
hold on;
pause(1)
O2A = 0.2; O4B = 0.7; BC = 0.6;
O4O2 = 0.3; O4G4 = 0.4; BG5 = 0.3; yC = 0.9;
m3= 0.5; m4= 6; m5= 4; m6= 2; J4=10; J5= 6; F=1000;
T=0; J2= 100; w= 75*2*pi/60; J3=0;
p0 = plot(0,0,'o'); hold on;
p2 = plot(0,O4O2,'o');
for n = 1:length(yi)
xA = O2A*cos(yi(1,n));
yA = O4O2+O2A*sin(yi(1,n));
xG4= O4G4*cos(yi(3,n));

yG4= O4G4*sin(yi(3,n));
xB = O4B*cos(yi(3,n));
yB = O4B*sin(yi(3,n));
xG5= O4B*cos(yi(3,n))+BG5*cos(yi(5,n));
yG5= O4B*sin(yi(3,n))+BG5*sin(yi(5,n));
xC = O4B*cos(yi(3,n))+BC*cos(yi(5,n));
yC = O4B*sin(yi(3,n))+BC*sin(yi(5,n));
p3 = plot([0 xA],[O4O2 yA],'b-');
p4 = plot([0 xB],[0 yB], 'b-');
p5 = plot([xB xC],[yB yC],'b-'); 
% pause time speed of simulation
pause(0.001)
delete(p3)
delete(p4)
delete(p5)
xlim([-1 0.6])
ylim([-0.5,1.5])
end
function [ydot, labda, Slider3] = calcydot(yi)
% yi, current value -> to ydot velocity and acc
g = yi(1);
gd = yi(2);
a = yi(3);
ad = yi(4);
b = yi(5);
bd = yi(6);
symb_Acc ;
gdd = Acc(1);
add = Acc(2);
bdd = Acc(3);
labda = Acc(4:5,1);
symb_slid;
Slider3 = slid;
%projection version with constrain
ydot = [gd; gdd; ad; add; bd; bdd];
end
%TMT method

function TMT_1()
syms g a b
syms gd ad bd gdd add bdd % for dx/dt I use xd etc
O2A = 0.2; O4B = 0.7; BC = 0.6;
O4O2 = 0.3; O4G4 = 0.4; BG5 = 0.3; yC = 0.9;
m3= 0.5; m4= 6; m5= 4; m6= 2; J4=10; J5= 6; F=1000;
T=0; J2= 100; w= 75*2*pi/60; J3=0;
%setting up kinematics.
q = [g;a;b];
qd = [gd;ad;bd];
qdd = [gdd; add; bdd];
xA = O2A*cos(g);
yA = O4O2+O2A*sin(g);
xG4= O4G4*cos(a);
yG4= O4G4*sin(a);
xG5= O4B*cos(a)+BG5*cos(b);
yG5= O4B*sin(a)+BG5*sin(b);
xC = O4B*cos(a)+BC*cos(b);
%making Ti,k matrix
Ki = [g;xA;yA;a;xG4;yG4;a;xG5;yG5;b;xC];%all kinematics combined
Ti = jacobian(Ki, q);
Tiv = Ti*qd; %velocity
gk = jacobian(Tiv,q)*qd;
Tacc = jacobian(Tiv,qd)*qdd + jacobian(Tiv,q)*qd;
%making Mass and gravity matrix and vector
M = diag([J2 m3 m3 J3 m4 m4 J4 m5 m5 J5 m6]);
% Constrains
Cy = yC - O4B*sin(a)-BC*sin(b);
ConA = O2A*cos(g)-sqrt((xA^2+yA^2))*cos(a);
Constra = [Cy;ConA];
CCd = jacobian(Constra,q);
CCdd = jacobian(CCd*qd,q)*qd;

%Forces
fi = zeros(11,1);
fi(1) = T;
fi(11) = F;
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
a_start = atan((O4O2 + sin(g)*O2A)/(cos(g)*O2A));
ad_start = jacobian(a_start,q)*qd;
b_start = pi-asin((yC-O4B*sin(a_start))/BC);
bd_start = jacobian(b_start,q)*qd;
angels_start = [a_start; ad_start; b_start; bd_start];
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
