%reaction rate parameters

Kr=100;

gamma_r=10;
kp=5;
%gamma_p=2;
%gamma_fb=6;
omega=100;
kr=Kr/omega;
N=1;
n1=1;

MAK=1;
X0=[0;0];

W_star=[];
C=[Kr,gamma_r,kp];
S=[1 -1 0; 0 0 1];
S_bis=[0 1 1; 0 0 0];
resample=1;
N_tot=10000;
delta=0.1;
t=0:delta:1;
tf=t(end);





x=zeros(2,length(t));
x(1,1)=10^(-3);


for k=1:length(t)-1
    x(1,k+1)=x(1,k)+(delta*(kr-(gamma_r*x(1,k))));
    x(2,k+1)=x(2,k)+(delta*(kp*x(1,k)));
end

%M_tot=x(1,:)+(sqrt(1/omega)*m);
    

N_state=500;

Hidden_Species=Hidden_State(0: N_state);
Shape_State_Space_Hidden=size(Hidden_Species);
State_Hidden=Shape_State_Space_Hidden(1);
[c,index]=intersect(Hidden_Species,X0(1:n1)','rows');
p0_1=zeros(State_Hidden,1);
p0_1(index)=1;

[t_1,t_ob,delta_1,match,X_tot_prev,X_1,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK);
Y_omega=Y./omega;
[V_tot_1,w_tot,V_jump,w_jump,match_1,match_2,resampling,pt,E_PF,Var_PF,SD_PF,E_pf] = particle_filter(t_ob,Y,p0_1,C,tf,N_tot,Hidden_Species,S,S_bis,resample,n1,W_star,MAK);
E_PF=E_PF./omega;
Var_PF=Var_PF./(omega)^2;
SD_PF=SD_PF./omega;
X_1=X_1./omega;
X_1_tot=X_tot_prev(n1,:)./omega;
%time flow of the process


[T,F,jump_times,E_FSP,Var_FSP,SD_FSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0_1,C,Hidden_Species,delta_1,S,S_bis,n1,W_star,MAK);
E_FSP=E_FSP./omega;
SD_FSP=SD_FSP./omega;
X_1=X_1./omega; %exact trajectory of the hidden process at the observation process jump times
X_1_tot=X_tot_prev(n1,:)./omega; %exact trajectory of the hidden process 

%%
shape_P=size(F);

if mod(length(t_ob),2) ~= 0
  dim=floor(length(t_ob)/2+1);
  medium=jump_times(2)-jump_times(1);
  if dim>2
  for j=2:1:dim
      medium=[medium jump_times(j+1)-jump_times(j)];
  end
  %medium=[medium jump_times(end)-jump_times(length(jump_times)-1)];
  end
else
    dim=length(t_ob)/2;
    if dim>2
    medium=jump_times(2)-jump_times(1);
    for j=2:1:dim-1
        medium=[medium jump_times(j+1)-jump_times(j)];
    end
    end
end
sampling=floor(mean(medium));
E_tot_bis=E_tot(1,:);
T_tot_bis=T(:,1);
SD_tot_bis=SD_tot(1,:);
count=1+sampling;
while count <= shape_P(2)
    E_tot_bis=[E_tot_bis E_tot(count,:)];
    T_tot_bis=[T_tot_bis T(:,count)];
    SD_tot_bis=[SD_tot_bis SD_tot(count,:)];
    count=count+sampling;
end
    


E_FFSP=E_tot_bis./omega;
SD_FFSP=SD_tot_bis./omega;
    


E=zeros(N,length(t));
Sigma=zeros(N,length(t));
SD_LNA=zeros(N,length(t));
b=1;
for j=1:2:(2*N)-1 %generation of N trajectories for each X component

%deterministic limit and fluctuation process initiliation 

mu=zeros(1,length(t));
sigma=zeros(1,length(t));



%Euler-Mayukama discretization in time for L1,L2,x,y,\mu,\Sigma
for i=1:length(t)-1
    %LNA approximation for X and Y
    %L(1,i+1)=L(1,i)+(delta*((-gamma_r))*L(1,i))+(sqrt(delta*kr)*normrnd(0,1))-(sqrt(delta*(gamma_r*x(1,i)))*normrnd(0,1));
    %L(2,i+1)=L(2,i)+(delta*kp*L(1,i))+(sqrt(delta*kp*x(1,i))*normrnd(0,1));
    
    %conditional mean evolution equation
    mu(i+1)=mu(i)+(delta*(-gamma_r))+(sigma(i)*(x(1,i))^(-1))*((Y_omega(i+1)-Y_omega(i))-(x(2,i+1)-x(2,i))-(kp*mu(i)*delta));
    %conditional variance evolution equation
    D=[sqrt(kr); -sqrt(gamma_r*x(1,i))];
    sigma(i+1)=sigma(i)+delta*((2*(-gamma_r)*sigma(i))+(transpose(D)*D)-((sigma(i)^(2))*kp*(x(1,i))^(-1)));
    
   

end
%LNA for X
%X(j,:)=x(1,:)+(L(1,:)/sqrt(omega));
%LNA for Y
%X(j+1,:)=x(2,:)+(L(2,:)/sqrt(omega));
%LNA for the filter E(t)=x(t)+sqrt(1/omega)*mu(t)
E(b,:)=x(1,:)+(mu/sqrt(omega));
Sigma(b,:)=sigma./omega;
SD_LNA=sqrt(Sigma);
b=b+1;


end


%%
f=figure;
f.Units='points';
f.OuterPosition=[10 10 1000 450];
freq=2;
subplot(1,2,1)
stairs(t_ob,Y_omega*100,'m','LineWidth',3)
xlabel('t (s)')
ylabel('Y(t)')
title('Observation Process')
xlim([0 t_ob(end)])
%xticks = get(gca,'ytick'); 
%scaling  = 10^2; 
%newlabels = arrayfun(@(y) sprintf('%.1f', scaling * y), yticks, 'un', 0);
Ax = gca;
Ax.XAxis.Exponent = 0;
Ax.YAxis.Exponent = 2;
set(gca,'FontSize',20)

subplot(1,2,2)
stairs(t_1,X_1_tot*100,'b','LineWidth',3)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_PF(([1 [1:freq:end]]),1)'*100,SD_PF(([1 [1:freq:end]]),1)*100,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',3)
hold on
errorbar(t,E*100,SD_LNA*100,'--','LineWidth',2)
hold on
%errorbar(T_tot_bis([1 [1:freq:end]]),E_FFSP([1 [1:freq:end]])*100,SD_FFSP([1 [1:freq:end]])*100,'r--','LineWidth',3)
errorbar(t_ob([1 [1:freq:end]]),E_FSP([1 [1:freq:end]])*100,SD_FSP([1 [1:freq:end]])*100,'r--','LineWidth',3)
% hold off
xlim([0 t_ob(end)])
%hold on
%stairs(t_ob,Y(1,:)./omega,'LineWidth',2)
xlabel('t (s)')
ylabel('Molecular Counts')
%legend('Exact Trajectory S(t) (Hidden Species)', 'Filtering Estimation of the S(t) (PF)' , 'Filtering Estimation of the S(t) (LNA-Kalman Filter)' ,'Observation Process A(t) trajectory' )
%legend('Exact', 'Particle Filter' , 'Kalman','FFSP','Location','northwest')
legend('Exact', 'Particle Filter' , 'Kalman','FFSP','Location','northwest')
title('Hidden Process')
%xticks = get(gca,'ytick'); 
%scaling  = 10^2; 
%newlabels = arrayfun(@(y) sprintf('%.1f', scaling * y), yticks, 'un', 0);
Ax = gca;
Ax.XAxis.Exponent = 0;
Ax.YAxis.Exponent = 2;
set(gca,'FontSize',20)


