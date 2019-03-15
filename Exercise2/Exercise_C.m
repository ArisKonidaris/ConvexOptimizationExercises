% Vissarion Konidaris
% Convex Optimization Course
% 1/4/2018

clear ; close all; clc

format long;
disp('Minimization the convex function');
disp(' c''x-sum(log(b-Ax)).');

in=0;
while(in~=1 && in~=2 && in~=3)
  disp('Give the dimensions of the problem.');
  disp('Indicative value pairs are (n,m) = (2,50), ');
  disp('(50,200),(300,800)\nwhere m the number of the ');
  disp('logarithms and n the dimension of input x.');
  disp('Type 1 for (2,50), 2 for (5,200), 3 for (300,800).');
  in=input('n : ');
  if(in==1)
    n=2;
    m=20;
  elseif(in==2)
    n=50;
    m=200;
  else
    n=300;
    m=800;
  end
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
A=2*limit*rand(m,n)-limit;
b=limit*rand(m,1); % A random vector b.
c=2*limit*rand(n,1)-limit; % A random vector c.

cvx_begin
    variable cv_x(n)
    minimize(c'*cv_x-sum(log(b-A*cv_x)))
cvx_end

optimal_value = c'*cv_x-sum(log(b-A*cv_x));

xaxis=linspace(cv_x(1,1)-1,cv_x(1,1)+1,250);
yaxis=linspace(cv_x(2,1)-1,cv_x(2,1)+1,250);
[X,Y]=meshgrid(xaxis,yaxis);
func = ones(size(X,1),size(X,2))*100;
if(n==2)
    for i=1:size(X,1)
        for j=1:size(X,2)
            log_argument = b-A*[X(i,j); Y(i,j)];
            feasible=1;
            for k=1:size(log_argument,1)
                if(log_argument(k,1)<=0)
                    feasible=0;
                    break;
                end
            end
            if(feasible==1)
                func(i,j) = c'*[X(i,j); Y(i,j)]-...
                            sum(log(log_argument));
            end
        end
    end
    
    figure(1);
    mesh(X,Y,func,'LineWidth',2);
    rotate3d on;
    axis tight;
    title('f:R^2+ -> R    f(x) = c^Tx-sum(log(b-Ax))');
    zlabel('f(x1,x2)');
    xlabel('x1');
    ylabel('x2');
    xaxis([]);
    grid on;

    figure(2);
    contourf(X,Y,func,10,'Linewidth',2);
    axis tight;
    colormap hot;
    title('f:R^2+ -> R    f(x) = c^Tx-sum(log(b-Ax))');
    xlabel('x1');
    ylabel('x2');
    grid on;

    figure(3);
    contour3(X,Y,func,'Linewidth',2);
    rotate3d on;
    axis tight;
    colormap hot;
    title('f:R^2+ -> R    f(x) = c^Tx-sum(log(b-Ax))');
    xlabel('x1');
    ylabel('x2');
    grid on;
end

% Gradient Descent using back-tracking line search.
disp('Gradient Descent using back-tracking line search.');

aplha=0.5;
beta=1.;
while(aplha>=0.5 || aplha<=0.)
  disp('Give a value for the parameter aplha.');
  disp('Acceptable values are in the open interval (0,0.5).');
  aplha=input('aplha : ');
  disp('');
end
while(beta>=1. || beta<=0.)
  disp('Give a value for the parameter beta.');
  disp('Acceptable values are in the open interval (0,1).');
  beta=input('beta : ');
  disp('');
end

xo =zeros(n,1);
x=xo;
k1=0;
grad=c+sum(A.*(b-A*x).^(-1))';
grad_norm=norm(grad);
traj1=[xo];
while(grad_norm>e)
  t=1;
  while(1)
      feasible=1;
      point = b-A*(x-t*grad);
      for i=1:size(point,1)
        if(point(i,1)<=0)
            feasible=0;
            break;
        end
      end
      if(feasible==0)
          disp(['Not in domain at loop ',num2str(k1+1)]);
          t=beta*t;
          continue;
      end
      break;
  end
  f=c'*x-sum(log(b-A*x));
  while(c'*(x-t*grad)-sum(log(b-A*(x-t*grad)))...
        > f-aplha*t*grad_norm^2)
    t=beta*t;
  end
  x=x-t*grad;
  k1=k1+1;
  grad=c+sum(A.*(b-A*x).^(-1))';
  grad_norm=norm(grad);
  traj1 = [traj1 x];
end
disp(['Actual iterations : ',num2str(k1)]);
disp(['Error : ',num2str(norm(cv_x-x))]);
disp(['Optimal value (cvx_optval): ',...
      num2str(optimal_value)]);
disp(['Optimal value (back-tracking): ',...
      num2str(c'*x-sum(log(b-A*x)))]);
disp('');
disp('');

% Newtons method
disp('Newtons method.');
x=xo;
k2=0;
grad=c+sum(A.*(b-A*x).^(-1))';
hessian=zeros(n,n);
for i=1:size(A,1)
    hessian=hessian+(A(i,:)'*A(i,:)).*...
            (norm(b(i,1)-A(i,:)*x).^2).^(-1);
end
traj2=[xo];
while(1)
  t=1;
  delta=-inv(hessian)*grad;
  lambda_squared=-grad'*delta;
  if(0.5*lambda_squared<=e)
      break;
  end
  while(1)
      feasible=1;
      point = b-A*(x+t*delta);
      for i=1:size(point,1)
        if(point(i,1)<=0)
            feasible=0;
            break;
        end
      end
      if(feasible==0)
          disp(['Not in domain at loop ',num2str(k2+1]));
          t=beta*t;
          continue;
      end
      break;
  end
  f=c'*x-sum(log(b-A*x));
  while(c'*(x+t*delta)-sum(log(b-A*(x+t*delta)))...
        > f-aplha*t*lambda_squared)
    t=beta*t;
  end
  x=x+t*delta;
  k2=k2+1;
  grad=c+sum(A.*(b-A*x).^(-1))';
  hessian=zeros(n,n);
  for i=1:size(A,1)
    hessian=hessian+(A(i,:)'*A(i,:)).*...
            (norm(b(i,1)-A(i,:)*x).^2).^(-1);
  end
  traj2 = [traj2 x];
end
disp(['Actual iterations : ',num2str(k2)]);
disp(['Error : ',num2str(norm(cv_x-x))]);
disp(['Optimal value (cvx_optval): ',...
      num2str(optimal_value)]);
disp(['Optimal value (Newtons): ',...
      num2str(c'*x-sum(log(b-A*x)))]);

v1=zeros(1,size(traj1,2));
for i=1:size(traj1,2)
    v1(1,i)=c'*traj1(:,i)-sum(log(b-A*traj1(:,i)));
end
iterations1 = linspace(1,k1+1,k1+1);

v2=zeros(1,size(traj2,2));
for i=1:size(traj2,2)
    v2(1,i)=c'*traj2(:,i)-sum(log(b-A*traj2(:,i)));
end
iterations2 = linspace(1,k2+1,k2+1);

figure(4);hold on;
semilogy(iterations1,v1-optimal_value,'linewidth',3);
semilogy(iterations2,v2-optimal_value,'linewidth',3);
title('Logarithmic Error');
xlabel('k');
ylabel('log(f(x_k)-p_*)');
legend('Descent with back-tracking line search','Newtons Method');
grid on;
