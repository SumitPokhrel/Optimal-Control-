
clear
clc
% Sumit Pokhrel  
% 03/20/2013
% Problem 2. a) 
% Optimal Control Problem using continuous-time necessary condition 
A = [0 1;0 0];
B = [0;1];
Q = [1 0;0 1];
Sf =[1 0;0 1];
R =1; 
H = [A -B*inv(R)*B'; -Q -A'];
He = expm (H);
Q11 = He(1:2,1:2);
Q12 = He(1:2,3:4);
Q21 = He(3:4,1:2);
Q22 = He(3:4,3:4);
alpha =1; % for counter

x0 = [1;1];
lamda0=inv(Q22-Sf*Q12)*(Sf*Q11-Q21)*x0;

x = x0; % initializing
lamda = lamda0;

  for t = 0:0.05:1
     x_and_lamda(alpha,:) = (expm (H*t)*[x;lamda])';
     u(alpha) = -R^(-1)*B'*x_and_lamda(alpha,3:4)';
     Hamiltonian(alpha)= 1/2*(x_and_lamda(alpha,1:2)*Q*x_and_lamda(alpha,1:2)'+u(alpha)'*R*u(alpha))+x_and_lamda(alpha,3:4)*(A*x_and_lamda(alpha,1:2)'+B*u(alpha));
     alpha = alpha+1;
  end 
  t=0:0.05:1; % since for loop cannot build a vector 

%Problem 2 b)
% Optimal Control Problem using sympletic-Euler method  

alpha =1; %for counter
x_and_lamda1(1,:) = [x;lamda];
for t = 0:0.05:1.0
    x_and_lamda1(alpha+1,1:2) = x_and_lamda1(alpha,1:2)'+(-B*R^(-1)*B'*x_and_lamda1(alpha,3:4)'+A*x_and_lamda1(alpha,1:2)')*0.05;
    x_and_lamda1(alpha+1,3:4) = x_and_lamda1(alpha,3:4)'+(-Q*x_and_lamda1(alpha,1:2)'- A'*x_and_lamda1(alpha,3:4)')*0.05;
    u1(alpha)  = -R^(-1)*B'*x_and_lamda1(alpha,3:4)';
    Hamiltonian1(alpha)= 1/2*(x_and_lamda1(alpha,1:2)*Q*x_and_lamda1(alpha,1:2)'+u1(alpha)'*R*u1(alpha))+x_and_lamda1(alpha,3:4)*(A*x_and_lamda1(alpha,1:2)'+B*u1(alpha));
    alpha = alpha+1;
end
 t=0:0.05:1.0; % since for loop cannot build a vector 
 
figure
plot(t,u,t,u1);
legend ('Continuous time','Symplectic Euler');
ylabel('u');
xlabel('time (sec)');
title('Control  vs. Time');

figure
plot(t,Hamiltonian,t,Hamiltonian1);
legend ('Continuous time','Symplectic Euler');
ylabel('H');
xlabel('time (sec)');
title('Hamiltonian vs. Time');









