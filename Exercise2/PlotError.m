function PlotError(iter,traj,P,q,x_star)
  
    tr_traj=traj';
    tr_traj=traj';
	v=zeros(1,size(traj,2));
	for i=1:size(traj,2)
		v(1,i)=0.5*tr_traj(i,:)*P*traj(:,i)+q'*traj(:,i);
	end
    p_star=0.5*x_star'*P*x_star+x_star'*q;
    func = log(v-p_star);
    iter = linspace(1,iter+1,iter+1);
  
    plot(iter,func,'linewidth',3);
    title('Logarithmic Error');
    xlabel('k');
    ylabel('log(f(x_k)-p_*)');
    grid on;
  
    return
end