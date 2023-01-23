clc
clear all
syms  m_1 m_2 l_2 l_1 g M 
%Case 1: When x of the system is only measurable or observable
A=[0 1 0 0 0 0;
   0 0 (-g*m_1)/M 0 (-g*m_2)/M 0;
   0 0 0 1 0 0;
   0 0 (-g*(M+m_1))/(l_1*M) 0 (-g*m_2)/(M*l_1) 0;
   0 0 0 0 0 1;
   0 0 (-g*m_1)/(M*l_2) 0 (-g*(M+m_2))/(M*l_2) 0];
C_1=[1 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];
O_1 = [C_1 ; C_1*A ; C_1*A^2 ; C_1*A^3 ; C_1*A^4 ; C_1*A^5];
lat_1 = O_1;
lat__1 = latex(lat_1);
disp("Observability of system when only x is measureable is given by Rank of the observability matrix");
disp("Rank using manual derivation of observability matrix is given by:")
disp(rank(O_1));
if(rank(O_1)<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");

%% Substituting values of M, m_1, m_2, l_1, l_2
%Case 1: When x of the system is only measurable or observable with values
%substituted
A_0 = double(subs(A,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
disp("Rank using Matlab function of observability is given as(after substituting variables: ");
Observ_using_matlab = rank(obsv(A_0,C_1));
disp(Observ_using_matlab);
O1_ = obsv(A_0,C_1);
lat_ = latex(sym(O1_));
if(Observ_using_matlab<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");
%% Case 2: When theta_1 and theta_2 of the system is only measurable or observable
C_2=[0 0 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 1 0];
O_2 = [C_2 ; C_2*A ; C_2*A^2 ; C_2*A^3 ; C_2*A^4 ; C_2*A^5];
disp("Observability of system when only theta_1 and theta_2 are observed is given by Rank of the observability matrix");
disp("Rank using manual derivation of observability matrix is given by:")
disp(rank(O_2));
lat_2 = O_2;
lat__2 = latex(simplify(lat_2));
if(rank(O_2)<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");
%% Case2 : When theta-1 and theta-2 are observed with values of Mass and lengths substituted
A_1 = double(subs(A,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
disp("Rank using Matlab function of observability is given as(after substituting variables: ");
Observ_using_matlab_1 = rank(obsv(A_1,C_2));
disp(Observ_using_matlab_1);
O2_ = obsv(A_1,C_2);
lat_1_ = latex(simplify(sym(O2_)));
if(Observ_using_matlab_1<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");
%% Case3: When x and theta_2 are taken as feedback for the system
C_3=[1 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 1 0];
O_3 = [C_3 ; C_3*A ; C_3*A^2 ; C_3*A^3 ; C_3*A^4 ; C_3*A^5];
disp("Observability of system when only x and theta_2 are observed is given by Rank of the observability matrix");
disp("Rank using manual derivation of observability matrix is given by:")
disp(rank(O_3));
lat_3 = latex(O_3);
if(rank(O_3)<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");
%%  Case3: When x and theta_2 are taken as feedback for the system with values of mass and lengths substituted
A_2 = double(subs(A,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
disp("Rank using Matlab function of observability is given as(after substituting variables: ");
Observ_using_matlab_2 = rank(obsv(A_2,C_3));
disp(Observ_using_matlab_2);
lat__3 = latex(sym(obsv(A_2,C_3)));
if(Observ_using_matlab_2<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");
%% Case4 : when x, theta_1, theta_2 are observable
C_4=[1 0 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 1 0];
O_4 = [C_4 ; C_4*A ; C_4*A^2 ; C_4*A^3 ; C_4*A^4 ; C_4*A^5];
disp("Observability of system when only x, theta_1 and theta_2 are observed is given by Rank of the observability matrix");
disp("Rank using manual derivation of observability matrix is given by:")
disp(rank(O_4));
lat_4 = latex(O_4);
if(rank(O_4)<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");
%% Case4 : when x, theta_1, theta_2 are observable
A_3 = double(subs(A,[M m_1 m_2 l_1 l_2 g], [1000 100 100 20 10 9.8]));
disp("Rank using Matlab function of observability is given as(after substituting variables: ");
Observ_using_matlab_2 = rank(obsv(A_3,C_4));
disp(Observ_using_matlab_2);
lat__4 = latex(sym(obsv(A_3,C_4)));
if(Observ_using_matlab_2<6)
    disp("System is not observable using the given observable parameters");
else
    disp("System is observable since matrix is of full rank");
end
disp("______________________________________________________________________");