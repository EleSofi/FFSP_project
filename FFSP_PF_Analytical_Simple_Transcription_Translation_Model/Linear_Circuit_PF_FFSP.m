X0=[0;5];
%Reaction Rates
C=[1;5;1];
%Stoichiometry Matrix
S=[1 0 0; 0 1 -1];
S_bis=[0 0 0; 1 0 1];
N_state=200;
n1=1;
n2=length(X0)-n1;
[State_Space] = Hidden_State(0:N_state);
Shape_State_Space=size(State_Space);
State=Shape_State_Space(1);
[c,index]=intersect(State_Space,X0(1:n1)','rows');
p0=zeros(length(State_Space),1);
p0(index)=1;
N_tot=10000;
MAK=1;
W_star=[];
resample=1;
tf=5;


[t,t_ob,delta,match,X_tot_prev,X,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK);

[V_tot,w_tot,V_jump,w_jump,match_1,match_2,resampling,pt,E_pf,Var_pf,SD_pf,E_pf_bis] = particle_filter(t_ob,Y,p0,C,tf,N_tot,State_Space,S,S_bis,resample,n1,W_star,MAK);


[T,F,jump_times,E_FSP,Var_FSP,SD_FSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);
 %exact trajectory of the hidden process at the observation process jump times
X_1_tot=X_tot_prev(n1,:); %exact trajectory of the hidden process 
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
E_tot_bis=E_tot(n1,:);
T_tot_bis=T(:,1);
SD_tot_bis=SD_tot(1,:);
count=1+sampling;
while count <= shape_P(2)
    E_tot_bis=[E_tot_bis E_tot(count,:)];
    T_tot_bis=[T_tot_bis T(:,count)];
    SD_tot_bis=[SD_tot_bis SD_tot(count,:)];
    count=count+sampling;
end
    


E_FFSP=E_tot_bis;
SD_FFSP=SD_tot_bis;

tm=3;
lambda=Poisson_Rate(tm,t_ob,Y(end,:),C(1));
u=@(x)(lambda.^(x)*exp(-lambda))./factorial(x);
t_indices=find(t_ob < tm);


%% plot

f=figure;
f.Units='points';
f.OuterPosition=[10 10 1000 950];
freq=2;
subplot(2,2,1)
stairs(t_ob,Y(1,:),'m','LineWidth',2)
xlim([0 t_ob(end)])
%xlabel("t (s) "+newline+"   ")
xlabel("t (s)")
ylabel('Y(t)')
title('Observation Process')
set(gca,'FontSize',20)  

end_state=22;            
subplot(2,2,2)
plot(State_Space(1:end_state),u(State_Space(1:end_state)),'b','LineWidth',5)
hold on
plot(State_Space(1:end_state),pt(1:end_state,t_indices(end)+1),'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold on
plot(State_Space(1:end_state),F(1:end_state,jump_times(t_indices(end)+1)),'r','LineWidth',3)
hold off
title('Conditional Probability')
xlabel('Hidden Process Copy Number')
xlim([0 end_state])
legend('Exact','PF','FFSP')
set(gca,'FontSize',20)


subplot(2,2,3)
stairs(t,X_1_tot,'b','LineWidth',2)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_FSP([1 [1:freq:end]],1)',SD_FSP([1 [1:freq:end]],1),'r--','LineWidth',3)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_pf([1 [1:freq:end]],1)',SD_pf([1 [1:freq:end]],1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
hold off
title('Hidden Process')
%xlabel("t (s) "+newline+"   ")
xlabel(" t(s) ")
ylabel('Molecular Counts')
xlim([0 t_ob(end)])
legend('Exact Trajectory', 'FFSP' , 'PF','Location','northwest')
set(gca,'FontSize',20)  
%set(gca, 'xtick', []);

% p = get(gca, 'Position');
% % Increase the height of the first and third subplots by 10%
% p_diff = p(4) * 0.1;
% % Increase the height of the subplot, but this will keep the
%         % bottom in the same place
% p(4) = p(4) + p_diff;
% % So also move the subplot down to decrease the gap to the next
% % one.
% p(2) = p(2) - p_diff;
% set(gca, 'Position', p);




Rel_1=(abs(pt(1:end_state,t_indices(end)+1)-u(State_Space(1:end_state))))./u(State_Space(1:end_state));
Rel_2=(abs(F(1:end_state,jump_times(t_indices(end)+1))-u(State_Space(1:end_state))))./u(State_Space(1:end_state));
subplot(2,2,4)
%plot(State_Space(1:end_state),Rel_1,'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
semilogy(State_Space(1:end_state),Rel_1,'Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold on
%plot(State_Space(1:end_state),Rel_2,'r','LineWidth',3)
semilogy(State_Space(1:end_state),Rel_2,'r','LineWidth',3)
hold off
title('Relative Error')
xlabel('Hidden Process Copy Number')
yticks([10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^(0)])
xlim([0 end_state])
legend('PF','FFSP')
set(gca,'FontSize',20)



%set(gca, 'xtick', []);
% p = get(gca, 'Position');
% % Increase the height of the first and third subplots by 10%
% p_diff = p(4) * 0.1;
% % Increase the height of the subplot, but this will keep the
%         % bottom in the same place
% p(4) = p(4) + p_diff;
% % So also move the subplot down to decrease the gap to the next
% % one.
% p(2) = p(2) - p_diff;
% set(gca, 'Position', p);
% 
% 

