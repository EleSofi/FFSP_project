X0=[1;0;0;0];
%Reaction Rates
C=[1,5,8,2,4,1];
%Stoichiometry Matrix
%S=[-1 1 0 0 0 0; 1 -1 0 0 0 0; 0 0 1 0 -1 0 ; 0 0 0 1 0 -1];
S=[1 -1 0 0 0 0; -1 1 0 0 0 0; 0 0 1 -1 0 0; -1 1 0 0 1 -1];
%S_bis=[1 0 0 0 0 0;0 1 1 0 0 0; 0 0 0 1 1 0; 0 0 0 0 0 1];
S_bis=[0 1 0 0 0 0; 1 0 1 0 0 0 ; 0 0 0 1 1 0; 1 0 0 0 0 1];
N_state=100;
n1=3;
A=zeros(N_state+1,n1);
A(:,n1)=0:N_state;
A(:,1)=1;
B=ones(N_state+1,n1);
B(:,1)=0;
B(:,n1)=0:N_state;
State_Space=vertcat(A,B);

n2=length(X0)-n1;

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
[T,F,jump_times,E_FSP,Var_FSP,SD_FSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);
[V_tot,w_tot,V_jump,w_jump,match_1,match_2,resampling,pt,E_pf,Var_pf,SD_pf,E_pf_bis] = particle_filter(t_ob,Y,p0,C,tf,N_tot,State_Space,S,S_bis,resample,n1,W_star,MAK);



X_1_tot=X_tot_prev(1:n1,:); %exact trajectory of the hidden process 
shape_P=size(F);

if mod(length(t_ob),2) ~= 0
  dim=floor(length(t_ob)/2+1);
  medium=jump_times(2)-jump_times(1);
  if dim>2
  for j=3:2:dim-2
      medium=[medium jump_times(j+1)-jump_times(j)];
  end
  medium=[medium jump_times(end)-jump_times(length(jump_times)-1)];
  end
else
    dim=length(t_ob)/2;
    if dim>2
    medium=jump_times(2)-jump_times(1);
    for j=3:2:dim
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
    E_tot_bis=[E_tot_bis; E_tot(count,:)];
    T_tot_bis=[T_tot_bis T(:,count)];
    SD_tot_bis=[SD_tot_bis ;SD_tot(count,:)];
    count=count+sampling;
end
    


E_FFSP=E_tot_bis;
SD_FFSP=SD_tot_bis;

%%
thresh=0.9;
for i=1:length(t_ob)
    if E_FSP(i,2)+SD_FSP(i,2) >= 1
        SD_FSP(i,2)=SD_FSP(i,2)-thresh;
    end
end



%% plot
f=figure;
f.Units='points';
f.OuterPosition=[10 10 2000 450];
freq=4;




subplot(1,3,1)
stairs(t_ob,Y(1,:),'m','LineWidth',2)
xlabel('t (s)')
ylabel('Y(t)')
xlim([0 t_ob(end)])
title('Observation Process')
set(gca,'FontSize',20)
%set(gca, 'xtick', []);


subplot(1,3,2)
stairs(t,X_1_tot(2,:),'b','LineWidth',2)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_FSP([1 [1:freq:end]],2)',SD_FSP([1 [1:freq:end]],2)','r--','LineWidth',3)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_pf([1 [1:freq:end]],2)',SD_pf([1 [1:freq:end]],1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold off
title('Hidden Process (Activated Gene)')
xlim([0 t_ob(end)])
ylim([0 1.5])
xlabel('t (s)')
ylabel('Molecular Counts')
legend('Exact Trajectory', 'FFSP' ,'PF','Location','northwest','Orientation','Horizontal')
set(gca,'FontSize',20)

% set(gca, 'xtick', []);
% 
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

subplot(1,3,3)
stairs(t,X_1_tot(3,:),'Color',[0.3010 0.7450 0.9330],'LineWidth',2)
hold on
%errorbar(T_tot_bis([1 [1:freq:end]]),E_FFSP([1 [1:freq:end]],3)',SD_FFSP([1 [1:freq:end]],3)','r--','LineWidth',3)
errorbar(t_ob([1 [1:freq:end]]),E_FSP([1 [1:freq:end]],3)',SD_FSP([1 [1:freq:end]],3)','r--','LineWidth',3)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_pf([1 [1:freq:end]],3)',SD_pf([1 [1:freq:end]],3),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold off
ylim([0 7])
xlim([0 t_ob(end)])
title('Hidden Process (mRNA)')
xlabel('t (s)')
ylabel('Molecular Counts')
legend('Exact Trajectory', 'FFSP' ,'PF','Location','northwest','Orientation','Horizontal')
set(gca,'FontSize',20) 
% %set(gca, 'xtick', []);
% p = get(gca, 'Position');
% % Increase the height of the first and third subplots by 10%
% p_diff = p(4) * 0.1;
% % Increase the height of the subplot, but this will keep the
% % bottom in the same place
% p(4) = p(4) + p_diff;
% % So also move the subplot down to decrease the gap to the next
% % one.
% p(2) = p(2) - p_diff;
% set(gca, 'Position', p);
