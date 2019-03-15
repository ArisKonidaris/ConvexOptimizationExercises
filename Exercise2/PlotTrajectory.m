function Plottraj(iter,traj,P,q,fig)
  
  tr_traj=traj';
  v=zeros(1,size(traj,2));
  for i=1:size(traj,2)
	v(1,i)=0.5*tr_traj(i,:)*P*traj(:,i)+q'*traj(:,i);
  end
  
  v = sort(v,'ascend');
  
  [m,idx]=max(abs(traj(1,iter+1)-traj(1,1:iter)));
  distx=abs(traj(1,iter+1)-traj(1,idx))+0.1;
  [m,idx]=max(abs(traj(2,iter+1)-traj(2,1:iter)));
  disty=abs(traj(2,iter+1)-traj(2,idx))+0.1;
  x=linspace(traj(1,iter+1)-distx,traj(1,iter+1)+distx,150);
  y=linspace(traj(2,iter)-disty,traj(2,iter)+disty,150);
  [X,Y]=meshgrid(x,y);
  func=0.5*(P(1,1).*(X.^2)+P(1,2).*X.*Y+P(2,1).*X.*Y+P(2,2).*(Y.^2))...
       +q(1,1).*X+q(2,1).*Y;
         
  figure(fig);
  contourf(X,Y,func,v);hold on;
  colormap hot;
  xlabel('x');
  ylabel('y');
  grid on;
  plot(traj(1,:),traj(2,:),'-d','MarkerSize',10,...
       'MarkerEdgeColor','b','MarkerFaceColor','b',...
       'linewidth',3);
  
  return
end