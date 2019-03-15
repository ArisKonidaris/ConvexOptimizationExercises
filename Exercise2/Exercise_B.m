clear ; close all; clc

format long;
disp('Minimization of a quadratic function');
disp(' with positive definite matrix.');
disp('');

% Dimensions of the problem. Indicative 
% values are n = 2, 50, 500, 1000.
n=0;
while(n~=2 && n~=50 && n~=500 && n~=1000)
  disp('Give the dimensions of the problem.');
  disp('Acceptable number of dimensions are ');
  disp('2, 50, 500 and 1000.');
  n=input('n : ');
  disp('');
end

% Condition Number
K=0;
while(K~=10 && K~=100 && K~=1000)
  disp('Give the condition number of matrix P.');
  disp('Acceptable number of condition numbers are ');
  disp('10, 100 and 1000.');
  K=input('K : ');
  disp('');
end

% The acceptable error of the optimizer.
e = -1.;
while(e<0.)
  disp('Give the acceptable error of the optimizer.');
  disp('Acceptable number of error are ');
  disp('the all the positive real numbers.');
  e=input('e : ');
  disp('');
end

% Number used for computing the range
% of the random values for the matrix A.
limit=100;

% The random matrix A constructed
% by n rand vectors of n dimensions.
A=2*limit*rand(n,n)-limit;

% Calculating the singular value decomposition of A.
[U,S,V]=svd(A);

% Constructing the eigenvalues for the matrix P.
lmin = randi(limit);
lmax = lmin*K;
L=diag([lmin;lmax;lmin+(lmax-lmin)*rand(n-2,1)]);

% Our positive definite matrix P with condition number K.
P=U*L*U';
q=2*limit*rand(n,1)-limit; % A random vector q.

% Closed for solution of the quadratic problem
sol=-inv(P)*q;

% Gradient Descent using exact line search.
disp('Gradient Descent using exact line search.');
xo = randn(n,1);
x = xo;
k1=0;
grad=P*x+q;
grad_norm=norm(grad);
traj1=x;
% Estimating the number of iterations through convergence analysis.
k1_e=K*log(((0.5*xo'*P*xo+xo'*q)-...
            (0.5*sol'*P*sol+sol'*q))...
            /e);
while(grad_norm>e)
  x = x-(grad_norm^2/(grad'*P*grad))*grad;
  k1=k1+1;
  grad=P*x+q;
  grad_norm=norm(grad);
  traj1 = [traj1 x];
end
disp(['Estimated iterations : ',num2str(k1_e)]);
disp(['Actual iterations : ',num2str(k1)]);
disp(['Error : ',num2str(norm(sol-x))]);
disp('')
disp('')


% Plot the trajectory of the minimisation algorithm.
if(n==2)
  PlotTrajectory(k1,traj1,P,q,1);
  title('Trajectory of descent using exact line search.');
end

% Gradient Descent using back-tracking line search.
disp('Gradient Descent using back-tracking line search.');

alpha=0.5;
beta=1.;
while(alpha>=0.5 || alpha<=0.)
  disp('Give a value for the parameter alpha.');
  disp('Acceptable values are in the open interval (0,0.5).');
  alpha=input('alpha : ');
  disp('');
end
while(beta>=1. || beta<=0.)
  disp('Give a value for the parameter beta.');
  disp('Acceptable values are in the open interval (0,1).');
  beta=input('beta : ');
  disp('');
end

x=xo;
k2=0;
grad=P*x+q;
grad_norm=norm(grad);
traj2=[x];
% Estimated number of iterations using convergence analysis.
k2_e=-(1/log(1-min(2*lmin*alpha,...
      (2*beta*alpha*lmin)/lmax)))*...
      log(((0.5*xo'*P*xo+xo'*q)-...
      (0.5*sol'*P*sol+sol'*q))...
      /e);
while(grad_norm>e)
  t=1;
  while(0.5*t^2*grad'*P*grad-t*grad_norm^2+0.5*x'*P*x+q'*x...
        > 0.5*x'*P*x+q'*x-alpha*t*grad_norm^2)
    t=beta*t;
  end
  x=x-t*grad;
  k2=k2+1;
  grad=P*x+q;
  grad_norm=norm(grad);
  traj2 = [traj2 x];
end
disp(['Estimated iterations : ',num2str(k2_e)]);
disp(['Actual iterations : ',num2str(k2)]);
disp(['Error : ',num2str(norm(sol-x))]);
disp('')
disp('')

% Plot the trajectory of the minimisation algorithm.
if(n==2)
  PlotTrajectory(k2,traj2,P,q,2);
  title('Trajectory of descent using back-tracking line search.');
end

% Plot log(f(x_k)-sol) against each iteration.
figure(3);hold on;
PlotError(k1,traj1,P,q,sol);
PlotError(k2,traj2,P,q,sol);
legend('exact line search','back-tracking line search');

% Plot theoretical and practical iterations
figure(4);
subplot(1,2,1);
hold on;
bar(categorical({'Estimated Iterations'}),k1_e,'FaceColor','blue');
bar(categorical({'Actual Iterations'}),k1,'FaceColor','red');
hold off;
title('Exact line search');
ax2=subplot(1,2,2);
hold on;
bar(categorical({'Estimated Iterations'}),k2_e,'FaceColor','blue');
bar(categorical({'Actual Iterations'}),k2,'FaceColor','red');
hold off;
title('Back-tracking line search');

disp('End of experiment.');
disp(['n : ',num2str(n)]);
disp(['K : ',num2str(K)]);
disp(['e : ',num2str(e)]);
disp(['alpha : ',num2str(alpha)]);
disp(['beta : ',num2str(beta)]);