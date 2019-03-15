% Vissarion Konidaris
% Convex Optimization Course
% 22/5/2018

clear ; close all; clc
format long;

disp('Minimization problem');
disp('dot(c,x)');
disp('subject to x>=0 and Ax=b');

in=0;
while(in~=1 && in~=2 && in~=3)
  disp('Give the dimensions of the problem.');
  disp('Indicative value pairs are (n,p) = (50,2), ');
  disp('(200,500),(800,300)\nwhere m the number of the ');
  disp('logarithms and n the dimension of input x.');
  disp('Type 1 for (50,2), 2 for (200,50), 3 for (800,300).');
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
  disp('Give the acceptable error of the optimizer.');
  disp('Acceptable number of error are ');
  disp('the all the positive real numbers.');
  e=input('e : ');
  disp('');
end

alpha=0.5;
beta=1.;
while(alpha>=0.5 || alpha<=0.)
  disp('Give a value for the parameter aplha.');
  disp('Acceptable values are in the open interval (0,0.5).');
  alpha=input('aplha : ');
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
c=rand(n,1);
b=A*x;


%%% Cvx solution %%%
%%%%%%%%%%%%%%%%%%%%

cvx_begin
    variable cv_x(n);
    minimize (c'*cv_x);
    subject to
        A*cv_x-b==0;
        -cv_x<=0;
cvx_end
cv_min_fun_val=c'*cv_x;


%%%%%%%%% Deasibility problem %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pseudoinverseA=pinv(A);
x0=pseudoinverseA*b;

s=max(-x0)+1;
mu=2;
tau=2;
feasible=0;
iter=0;
A_=[A zeros(size(A,1),1)];
while(1)
    while(1)
        temp = -(x0+s).^-1;
        grad_phi = [temp; sum(temp)];
        clear temp;
        grad = [zeros(n,1);  tau] + grad_phi;
        temp = (x0+s).^-2;
        hessian_phi = diag( [temp; sum(temp)] );
        hessian_phi(n+1,1:n) = temp';
        hessian_phi(1:n,n+1) = temp;
        clear temp;
        hessian_inv=pinv(hessian_phi);
        clear clear temp;
        
        v=-inv(A_*hessian_inv*A_')*A_*hessian_inv*grad;
        delta=-hessian_inv*(grad+A_'*v);
        lambda_squared=-grad'*delta;

        if(lambda_squared/2<=e)
            break;
        end

        t=1;
        while(1)
            in_domain=1;
            point = [x0; s] + t*delta;
            for i=1:n
                if(point(i,1)+point(n+1,1)<=0)
                    in_domain=0;
                    clear point;
                    break;
                end
            end
            if(in_domain==0)
                disp(['Not in domain at loop ',num2str(iter+1)]);
                t=beta*t;
                clear point;
                continue;
            end
            break;
        end
        f=tau*s-sum(log(x0+s));
        while( tau*(s+t*delta(n+1,1))-sum( log( (x0+t*delta(1:n,1))+...
               (s+t*delta(n+1,1)) ) )> f-alpha*t*lambda_squared)
            t=beta*t;
        end

        iter=iter+1;
        x0=x0+t*delta(1:n,1);
        s=s+t*delta(n+1,1);
        
        % Stopping earlier.
        if(s<0)
            feasible=1;
            break;
        end 
    end
    
    if((mu/tau)<e || feasible==1)
        break;
    end
    tau=tau*mu;
end



%%%%%% Barrier method starting from feasible point %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 2;
tau = 2;
iter = 0;
x1 = x0;
traj1 = [x1];
while(1)
    while(1)
        grad = tau*c-x1.^-1;
        hessian_inv=inv(diag(x1.^-2));
        v=-inv(A*hessian_inv*A')*A*hessian_inv*grad;
        delta=-hessian_inv*(grad+A'*v);
        lambda_squared=-grad'*delta;

        %disp(norm(cv_x-x))
        
        if(lambda_squared/2<=e)
            break;
        end

        t=1;
        while(1)
            in_domain=1;
            point = x1+t*delta;
            for i=1:size(point,1)
                if(point(i,1)<=0)
                    feasible=0;
                    clear point;
                    break;
                end
            end
            if(in_domain==0)
                disp(['Not in domain at loop ',num2str(iter+1)]);
                t=beta*t;
                clear point;
                continue;
            end
            break;
        end
        f=tau*c'*x1-sum(log(x1));
        while( tau*c'*(x1+t*delta)-sum(log(x1+t*delta) )...
            > f-alpha*t*lambda_squared)
            t=beta*t;
        end

        iter=iter+1;
        x1=x1+t*delta;
        traj1 = [traj1 x1];
    end
    
    if((mu/tau)<e)
        break;
    end
    tau=tau*mu;
end

%%%%%% Barrier method starting from infeasible point %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu=20;
tau=2;
epsilon_feas=0.001;
x2=x0;
l=rand(n,1);
v=rand(p,1);
en = x2'*l;
rho = [c-1*diag(ones(n,1))*l+A'*v;...
       diag(l)*x2-ones(n,1)*tau^-1;...
       A*x2-b];
norm_rho =norm(rho);
tau=mu*p/en;
traj2=[x2];
while(norm_rho>=epsilon_feas && en>=e)
    
    mat=zeros(2*n+p,2*n+p);
    mat(1:n,n+1:2*n) = -1*diag(ones(n,1));
    mat(1:n, 2*n+1:2*n+p) = A';
    mat(n+1:2*n,1:n) = diag(l)*diag(ones(n,1));
    mat(n+1:2*n,n+1:2*n) = diag(x2);
    mat(2*n+1:2*n+p,1:n) = A;
    
    deltas = -inv(mat)*rho;
     
    % Determinig step s.
    s=1;
    for i=n+1:2*n
        if(deltas(i,1)<0 && -l(i-n,1)/deltas(i,1)<s)
            s=-l(i-n,1)/deltas(i,1);
        end
    end
    
    % Backtracking s until f(x+)<0.
    s=0.99*s;
    while(1)
        feas=1;
        point = x2+s*deltas(1:n,1);
        for i=1:n
            if( -point(i,1) >=0 )
                feas=0;
                break;
            end
        end
        if(feas==1)
            break;
        end
        s=beta*s;
    end
    
    rtd_plus=c-1*diag(ones(n,1))*(l+s*deltas(n+1:2*n,1))+...
             A'*(v+s*deltas(2*n+1:2*n+p,1));
    rtc_plus=diag((l+s*deltas(n+1:2*n,1)))*(x2+s*deltas(1:n,1))-...
             ones(n,1)*tau^-1;
    rtp_plus=A*(x2+s*deltas(1:n,1))-b;
    rho_plus=[rtd_plus;rtc_plus;rtp_plus];
    norm_rho_plus=norm(rho_plus);
    while(norm_rho_plus^2 >(1-alpha*s)*norm_rho^2)
      s=beta*s;
      rtd_plus=c-1*diag(ones(n,1))*(l+s*deltas(n+1:2*n,1))+...
               A'*(v+s*deltas(2*n+1:2*n+p,1));
      rtc_plus=diag(l+s*deltas(n+1:2*n,1))*(x2+s*deltas(1:n,1))-...
               ones(n,1)*tau^-1;
      rtp_plus=A*(x2+s*deltas(1:n,1))-b;
      rho_plus=[rtd_plus;rtc_plus;rtp_plus];
      norm_rho_plus=norm(rho_plus);
    end
    x2 = x2+s*deltas(1:n,1);
    l = l+s*deltas(n+1:2*n,1);
    v = v+s*deltas(2*n+1:2*n+p,1);
    
    en = x2'*l;
    tau=mu*p/en;
    rho = [c-1*diag(ones(n,1))*l+A'*v;...
           diag(l)*x2-ones(n,1)*tau^-1;...
           A*x2-b];
    norm_rho = norm(rho);
    traj2 = [traj2 x2];
    
end

v1=zeros(1,size(traj1,2));
for i=1:size(traj1,2)
    v1(1,i)=c'*traj1(:,i);
end
iterations1 = linspace(1,size(traj1,2),size(traj1,2));

v2=zeros(1,size(traj2,2));
for i=1:size(traj2,2)
    v2(1,i)=c'*traj2(:,i);
end
iterations2 = linspace(1,size(traj2,2),size(traj2,2));

figure(1);hold on;
semilogy(iterations1,v1-v1(1,size(v1,2)),'linewidth',3);
semilogy(iterations2,v2-v2(1,size(v2,2)),'linewidth',3);
title('Logarithmic Error');
xlabel('k');
ylabel('log(f(x_k)-p_*)');
legend('Barrier Meth from feasible point','Primal-Dual Interior Point meth');
grid on;

disp(['Optimal value (cvx_optval): ',...
      num2str(c'*cv_x)]);
disp(['Optimal value (Barrier Mathod starting from feasible point): ',...
      num2str(c'*x1)]);
disp(['Optimal value (Primal-Dual Interior Point Method): ',...
      num2str(c'*x2)]);











