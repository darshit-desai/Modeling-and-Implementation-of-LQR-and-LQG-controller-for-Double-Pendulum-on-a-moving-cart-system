clc
clear all
syms  m_1 m_2 l_2 l_1 g M 
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
Controllability_Mat = [B, A*B, A^2*B, A^3*B, A^4*B, A^5*B];

r_c_mat=rank(Controllability_Mat);
disp("The rank of the matrix of controllability is = ");
r_c_mat

disp("The determinant of controllability matrix is=");
det(Controllability_Mat)


