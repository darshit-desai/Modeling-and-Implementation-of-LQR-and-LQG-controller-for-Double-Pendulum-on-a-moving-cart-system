clc
clear all
syms  m_1 m_2 l_2 l_1 g M 
%Initial Conditions
X0 = [1; 0; deg2rad(0); 0; deg2rad(0); 0]
A=[0 1 0 0 0 0;
   0 0 (-g*m_1)/M 0 (-g*m_2)/M 0;
   0 0 0 1 0 0;
   0 0 (-g*(M+m_1))/(l_1*M) 0 (-g*m_2)/(M*l_1) 0;
   0 0 0 0 0 1;
   0 0 (-g*m_1)/(M*l_2) 0 (-g*(M+m_2))/(M*l_2) 0];
B=[0;
    1/M;
    0;
    1/(M*l_1);
    0;
    1/(M*l_2)];
D=0;
C_1=[1 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];

Q = [1000 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];
R= 0.01;
%Substituting values of m1 l1 l2 and M in A
A1=double(subs(A,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
B1 = double(subs(B,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));

[K1,S,P] = lqr(A1, B1, Q, R);

sys = ss(A1-B1*K1,B1,C_1,D);

%Kalman Estimator Design
Bd = 0.01*eye(6); %disturbance
Vn = 0.001;     %Gaussian White Noise
[L,P,E] = lqe(A1,Bd,C_1,Bd,Vn*eye(3)); %Considering vector output: x(t)
Ac1 = A1-(L*C_1);
Xf = [20;0;0;0;0;0]
e_sys1 = ss(Ac1,[B1 L],C_1,0);
ts = 0:0.01:100;
[t,state1] = ode45(@(t,state)nonLinear1_LQG(t,state,-K1*(state-Xf),L),ts,X0);
figure();
hold on
plot(t,state1(:,1))
plot(t,state1(:,3))
plot(t,state1(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear LQG for output vector: x(t)')
legend('x','theta_1','theta_2')
hold off




%% LQG Non Linear
function ds1 = nonLinear1_LQG(t,x,f,l1)
    m1 = 100; m2 = 100; M=1000; L1 = 20; L2 = 10; g = 9.81;
    X_ = x(1);
    X_dot = x(2);
    theta_1 = x(3);
    theta_dot1 = x(4);
    theta_2 = x(5);
    theta_dot2 = x(6);
    ds1 = zeros(6,1);
    y1 = [X_;0;0];
    c_1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
    iter = l1*(y1-c_1*x);
    ds1(1) = X_dot+iter(1);
    ds1(2) = (f - m1*L1*sin(theta_1)*theta_dot1^2 - (m2*L2*sin(theta_2)*theta_dot2^2 - m1*g*sin(theta_1)*cos(theta_1) - m2*g*sin(theta_2)*cos(theta_2))/(M + m1 + m2 - m1*cos(theta_1)^2 - m2*cos(theta_2)^2))+iter(2);
    ds1(3) = theta_dot1+iter(3);
    ds1(4) = cos(theta_1)*ds1(2)/L1 - (g*sin(theta_1)/L1)+iter(4);
    ds1(5) = theta_dot2+iter(5);
    ds1(6) = cos(theta_2)*ds1(2)/L2 - (g*sin(theta_2)/L2)+iter(6);
end