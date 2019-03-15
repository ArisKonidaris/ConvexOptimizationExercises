#Vissarion Konidaris
#Convex Optimization Course
# 9/3/2018

clc;
clear;
clc;

%%%%%%%%%%%%%%%%%%% First Exercise %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start = 0.;
finish = 25.;
precision = 1000000;
domain = linspace(start,finish,precision);
func1 = 1./(1.+domain);

[FirstTaylor,SecondTaylor] = F1_Taylor_Appr(2.5,domain);
% Plot the function along with its first and
% second taylor approximations at point 2.5
figure(1);
pl = plot(domain,func1,'linewidth',2,'b',
          domain,FirstTaylor,'linewidth',2,'r',
          domain,SecondTaylor,'linewidth',2,'g');
xlabel('Domain');
ylabel('f(Domain)');
title('f:R+ -> R    f(x) = 1/(1+x)');
legend('f(5)','First Taylor Approximation of f at 2.5',
       'Second Taylor Approximation of f at 2.5');
grid on

[FirstTaylor,SecondTaylor] = F1_Taylor_Appr(5,domain);
% Plot the function along with its first and
% second taylor approximations at point 5
figure(2);
pl = plot(domain,func1,'linewidth',2,'b',
          domain,FirstTaylor,'linewidth',2,'r',
          domain,SecondTaylor,'linewidth',2,'g');
xlabel('Domain');
ylabel('f(Domain)');
title('f:R+ -> R    f(x) = 1/(1+x)');
legend('f(5)','First Taylor Approximation of f at 5',
       'Second Taylor Approximation of f at 5');
grid on

[FirstTaylor,SecondTaylor] = F1_Taylor_Appr(7.5,domain);
% Plot the function along with its first and
% second taylor approximations at point 7.5
figure(3);
pl = plot(domain,func1,'linewidth',2,'b',
          domain,FirstTaylor,'linewidth',2,'r',
          domain,SecondTaylor,'linewidth',2,'g');
xlabel('Domain');
ylabel('f(Domain)');
title('f:R+ -> R    f(x) = 1/(1+x)');
legend('f(5)','First Taylor Approximation of f at 7.5',
       'Second Taylor Approximation of f at 7.5');
grid on

[FirstTaylor,SecondTaylor] = F1_Taylor_Appr(10,domain);
% Plot the function along with its first and
% second taylor approximations at point 10
figure(4);
pl = plot(domain,func1,'linewidth',2,'b',
          domain,FirstTaylor,'linewidth',2,'r',
          domain,SecondTaylor,'linewidth',2,'g');
xlabel('Domain');
ylabel('f(Domain)');
title('f:R+ -> R    f(x) = 1/(1+x)');
legend('f(5)','First Taylor Approximation of f at 10',
       'Second Taylor Approximation of f at 10');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Second Exercise %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start = 0.;
finish = 20.;
precision = 250;
x1 = linspace(start,finish,precision);
x2 = x1;
[X,Y] = meshgrid(x1,x2);
func2 = 1./(1.+X.+Y);

% Plot the 3d plot of the second function using mesh
figure(5);
mesh(func2,'LineWidth',2);
rotate3d on;
axis([0 100 0 100 0 1]);
title('f:R^2+ -> R    f(x) = 1/(1+x1+x2)');
zlabel('f(x1,x2)');
xlabel('x1');
ylabel('x2');
grid on

% Plot the level sets of the second function using contour
figure(6);
contourf(func2,'LineWidth',2);
colormap hot
axis([0 100 0 100 0 1]);
title('f:R^2+ -> R    f(x) = 1/(1+x1+x2)    Level Sets');
xlabel('x1');
ylabel('x2');
grid on

[FirstTaylor,SecondTaylor] = F2_Taylor_Appr(2,2,X,Y);
% Common 3d plot of the second function with its
% first Taylor Approximation at point (2,2)
figure(7);hold on
mesh(func2,'LineWidth',2);
surf(FirstTaylor,'FaceColor','red','EdgeColor','none');
rotate3d on;
axis([0 150 0 150 0 1]);
title('f:R^2+ -> R    f(x) = 1/(1+x1+x2)');
zlabel('f(x1,x2)');
xlabel('x1');
ylabel('x2');
grid on

% Common 3d plot of the second function with its
% second Taylor Approximation at point (2,2)
figure(8);hold on
mesh(func2,'LineWidth',2);
surf(SecondTaylor,'FaceColor','red','EdgeColor','none');
rotate3d on;
axis([0 150 0 150 0 1]);
title('f:R^2+ -> R    f(x) = 1/(1+x1+x2)');
zlabel('f(x1,x2)');
xlabel('x1');
ylabel('x2');
grid on

[FirstTaylor,SecondTaylor] = F2_Taylor_Appr(5,5,X,Y);
% Common 3d plot of the second function with its
% first Taylor Approximation at point (5,5)
figure(9);hold on
mesh(func2,'LineWidth',2);
surf(FirstTaylor,'FaceColor','red','EdgeColor','none');
rotate3d on;
axis([0 150 0 150 0 1]);
title('f:R^2+ -> R    f(x) = 1/(1+x1+x2)');
zlabel('f(x1,x2)');
xlabel('x1');
ylabel('x2');
grid on

% Common 3d plot of the second function with its
% second Taylor Approximation at point (5,5)
figure(10);hold on
mesh(func2,'LineWidth',2);
surf(SecondTaylor,'FaceColor','red','EdgeColor','none');
rotate3d on;
axis([0 150 0 150 0 1]);
title('f:R^2+ -> R    f(x) = 1/(1+x1+x2)');
zlabel('f(x1,x2)');
xlabel('x1');
ylabel('x2');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% Fifth Exercise %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start = 0.01;
finish = 25.;
precision = 1000000;
domain = linspace(start,finish,precision);
func1 = domain.^3;
func2 = domain.^0.25;
func3 = domain.^(-2);

% Plot the function along with its first and
% second taylor approximations at point 2.5
figure(11);
pl = plot(domain,func1,'linewidth',2,'b',
          domain,func2,'linewidth',2,'r',
          domain,func3,'linewidth',2,'g');
xlabel('Domain');
ylabel('f(Domain)');
title('f:R++ -> R    f(x) = x^a');
legend('a=3','a=0.25','a=-2');
axis([-0.5 10 -2.5 25]);
grid on

start = -20.;
finish = 20.;
precision = 500;
x1 = linspace(start,finish,precision);
x2 = x1;
[X,Y] = meshgrid(x1,x2);
norm2 = (X.^2.+Y.^2).^(1/2);

% Plot the 3d plot of the second function using mesh
figure(12);
mesh(norm2,'LineWidth',2);
rotate3d on;
title('f:R^2 -> R    f(x) = ||X||');
zlabel('f(x1,x2)');
xlabel('x1');
ylabel('x2');
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%% Exercise 7 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = 200*rand(3,2)-100;
x = 200*rand(2,1)-100;
b = A*x;

precision = 15;
x1 = linspace(x(1,1)-abs(x(1,1))/1000,x(1,1)+abs(x(1,1))/1000,precision);
x2 = linspace(x(2,1)-abs(x(2,1))/1000,x(2,1)+abs(x(2,1))/1000,precision);
[X,Y] = meshgrid(x1,x2);

LeastSquares = (A(1,1)*X+A(1,2)*Y-b(1,1)).^2 ...
              +(A(2,1)*X+A(2,2)*Y-b(2,1)).^2 ...
              +(A(3,1)*X+A(3,2)*Y-b(3,1)).^2;
               
% Plot the 3d plot of the Least Squares Equation
figure(1);
mesh(X,Y,LeastSquares,'LineWidth',2);
rotate3d on;
title('f:R^2 -> R    f(x) = ||Ax-b||^2');
zlabel('f(x)');
xlabel('x1');
ylabel('x2');
grid on

% Plot the level sets of the Least Squares Equation 
figure(2);
contourf(X,Y,LeastSquares,'LineWidth',2);
colormap hot
title('f:R^2 -> R    f(x) = ||Ax-b||^2    Level Sets');
xlabel('x1');
ylabel('x2');
grid on

b+=randn(3,1); % Adding Gaussian noise to the b vector
LeastSquaresE = (A(1,1)*X+A(1,2)*Y-b(1,1)).^2 ...
               +(A(2,1)*X+A(2,2)*Y-b(2,1)).^2 ...
               +(A(3,1)*X+A(3,2)*Y-b(3,1)).^2;

% Plot the 3d plot of the Least Squares Equation
figure(3);
mesh(X,Y,LeastSquaresE,'LineWidth',2);
rotate3d on;
title('f:R^2 -> R    f(x) = ||Ax-be||^2');
zlabel('f(x)');
xlabel('x1');
ylabel('x2');
%xlim([X(1,1)-10 X(1,size(X,2))+10]);
%ylim([Y(1,1)-10 Y(size(Y,1),1)+10]);

grid on

% Plot the level sets of the Least Squares Equation 
figure(4);
contourf(X,Y,LeastSquaresE,'LineWidth',2);
colormap hot
title('f:R^2 -> R    f(x) = ||Ax-be||^2    Level Sets');
xlabel('x1');
ylabel('x2');
grid on

disp('x')
disp(x)
disp('X')
disp(X)
disp('Y')
disp(Y)
disp('LeastSquares')
disp(LeastSquares)
disp('LeastSquaresE')
disp(LeastSquaresE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










