% Vissarion Konidaris
% Convex Optimization Course
% 6/5/2018

clear ; close all; clc
format long;

disp('Minimization problem');
disp('-sum(log(x))');
disp('subject to Ax=b');

in=0;
while(in~=1 && in~=2 && in~=3)
  disp([newline 'Give the dimensions of the problem.']);
  disp('Indicative value pairs are (n,m) = (2,50), ');
  disp('(200,50),(800,300) where m the number of the ');
  disp('logarithms and n the dimension of input x.');
  disp('Type 1 for (2,50), 2 for (5,200), 3 for (800,300).');
  in=input('n : ');
  if(in==1)
    n=20;
    p=2;
  elseif(in==2)
    n=200;
    p=50;
  else
    n=800;
    p=300;
  end
  disp('');
end

% The acceptable error of the optimizer.
e = -1.;
while(e<0.)
  disp([newline 'Give the acceptable error of the optimizer.']);
  disp('Acceptable number of error are ');
  disp('the all the positive real numbers.');
  e=input('e : ');
  disp('');
end

aplha=0.5;
beta=1.;
while(aplha>=0.5 || aplha<=0.)
  disp([newline 'Give a value for the parameter aplha.']);
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

A=rand(p,n);
x=rand(n,1);
b=A*x;

%%% Cvx solution %%%

cvx_begin
    variable cv_x(n);
    minimize (-sum(log(cv_x)));
    subject to
        A*cv_x-b==0;
        -cv_x<0;
cvx_end
cv_min_fun_val=-sum(log(cv_x));
disp(cv_min_fun_val);


%%% Newtons' method starting from a feasible point %%%

% Gradient Descent using back-tracking line search.
disp([newline 'Newtons method starting from a feasible point.']);

% Finding a feasible starting point x0 via cvx
cvx_begin
    variable x0(n);
    minimize(0);
    subject to
        A*x0-b==0;
        -x0<0;
cvx_end

x=x0;
iter=0;
while(1)
    grad=-x.^-1;
    hessian=diag(x.^-2);
    hessian_inv=inv(hessian);
    w=-inv(A*hessian_inv*A')*A*hessian_inv*grad;
    delta=-hessian_inv*(grad+A'*w);
    lambda_squared=-grad'*delta;
   
    if(lambda_squared/2<=e)
        break;
    end
    
    t=1;
    while(1)
        feasible=1;
        point = x+t*delta;
        for i=1:size(point,1)
            if(point(i,1)<=0)
                feasible=0;
                break;
            end
        end
        if(feasible==0)
            disp(['Not in domain at loop ',num2str(iter+1)]);
            t=beta*t;
            continue;
        end
        break;
    end
    f=-sum(log(x));
    while(-sum(log(x+t*delta))...
        > f-aplha*t*lambda_squared)
        t=beta*t;
    end
    
    iter=iter+1;
    x=x+t*delta;
end
disp(['Actual iterations : ',num2str(iter)]);
disp(['Error : ',num2str(norm(cv_x-x))]);
disp(['Optimal value (cvx_optval): ',...
      num2str(cv_min_fun_val)]);
disp(['Optimal value (Newtons): ',...
      num2str(-sum(log(x)))]);
disp('')

%%% Newtons' method starting from a non feasible point %%%
disp([newline 'Primal-Dual algorithm starting from  a non feasible point.']);

x=ones(n,1);
v=rand(p,1);
iter=0;
stopping_norm=1.;
while(stopping_norm>e)
    grad=-x.^-1;
    hessian=diag(x.^-2);
    hessian_inv=inv(hessian);
    v_plus=-inv(A*hessian_inv*A')*(A*hessian_inv*grad-A*x+b);
    v_d=v_plus-v;
    x_d=-hessian_inv*(grad+A'*v_plus);
    
    t=1;
    while(1)
        feasible=1;
        point = x+t*x_d;
        for i=1:size(point,1)
            if(point(i,1)<=0)
                feasible=0;
                break;
            end
        end
        if(feasible==0)
            disp(['Not in domain at loop ',num2str(iter+1)]);
            t=beta*t;
            continue;
        end
        break;
    end
    f=sqrt(norm(grad+A'*v,2)^2+norm(A*x-b,2)^2);
    while(sqrt(...
          norm(-(x+t*x_d).^(-1)+A'*(v+t*v_d),2)^2+...
          norm(A*(x+t*x_d)-b,2)^2)...
        > (1-aplha*t)*f)
        t=beta*t;
    end
    
    iter=iter+1;
    x=x+t*x_d;
    v=v+t*v_d;
    stopping_norm=sqrt(norm(-x.^-1+A'*v,2)^2+norm(A*x-b,2)^2);
end
disp(['Actual iterations : ',num2str(iter)]);
disp(['Error : ',num2str(norm(cv_x-x))]);
disp(['Optimal value (cvx_optval): ',...
      num2str(cv_min_fun_val)]);
disp(['Optimal value (Newtons): ',...
      num2str(-sum(log(x)))]);
disp('');
  
%%% Cvx dual solution %%%
disp([newline 'Solving the Dual problem via cvx.']);
cvx_begin
    variable lag_mul(p);
    minimize (-sum(log(A'*lag_mul))-n+b'*lag_mul);
cvx_end
cvx_dual_sol=(A'*lag_mul).^-1;
cv_min_fun_val_dual=-sum(log(cvx_dual_sol));
disp(cv_min_fun_val_dual);
    