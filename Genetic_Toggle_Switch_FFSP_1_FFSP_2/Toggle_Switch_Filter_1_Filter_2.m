X0=[0;0];

%Reaction Rates
%C(1)=\alpha_1 C(2)=\alpha_2 C(3)=\beta C(4)=\gamma
%C=[50;16;2.5;1];
C=[16;14;1;1];
%C=[];
% alpha_1=5; these parameters work
% alpha_2=5;
% beta=1;
% gamma=1;
alpha_1=C(1);
alpha_2=C(2);
beta=C(3);
gamma=C(4);
%Stoichiometry Matrix
S=[1 -1 0 0; 0 0 1 -1];
S_bis=[0 1 0 0; 0 0 0 1];



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
MAK=0;
W_star{1}=@(X,Y)(alpha_1/(1+(Y).^(beta)));
W_star{2}=@(X,Y)X(1);
W_star{3}=@(X,Y)(alpha_2/(1+(X).^(gamma)));
W_star{4}=@(X,Y)Y;
resample=1;
tf=5;

[t,t_ob,delta,match,X_tot_prev,X,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK);
% stairs(t,X_tot_prev(2,:),'b')
% hold on
% stairs(t_ob,Y,'r')
%[T,F,jump_times,E_FSP,Var_FSP,SD_FSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);
%[V_tot,w_tot,V_jump,w_jump,match_1,match_2,resampling,pt,E_pf,Var_pf,SD_pf,E_pf_bis] = particle_filter(t_ob,Y,p0,C,tf,N_tot,State_Space,S,S_bis,resample,n1,W_star,MAK);

%FFSP1
[T,F,jump_times,E,Var,SD,Err_jump,E_tot,SD_tot,rho,eps] = FFSP_2_error(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);
%FFSP2
[T_1,P,jump_times_1,E_1,Var_1,SD_1,Err_jump_1,E_tot_1,SD_tot_1,Err_tot_1] = FFSP_2_error_new(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);



X_1_tot=X_tot_prev(n1,:); %exact trajectory of the hidden process 

%%
shape_P=size(F);

if mod(length(t_ob),2) ~= 0
  %dim=floor((length(t_ob)/2)+1);
  dim=length(t_ob)-1;
  medium=jump_times(2)-jump_times(1);
  if dim>2
  %for j=3:1:dim-2
  for j=2:1:dim
      medium=[medium jump_times(j+1)-jump_times(j)];
  end
  %medium=[medium jump_times(end)-jump_times(length(jump_times)-1)];
  end
else
    %dim=length(t_ob)/2;
    dim=length(t_ob);
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


E_tot_bis_1=E_tot_1(1,:);
T_tot_bis_1=T_1(:,1);
SD_tot_bis_1=SD_tot_1(1,:);
count_1=2;
%i=1+sampling;
i=2;
j=2;
tau=0.01;
while  (i <= length(T_1)) && (j <= length(T_tot_bis))
    if abs(T_1(i)-T_tot_bis(j)) <= tau
    E_tot_bis_1=[E_tot_bis_1; E_tot_1(i,:)];
    T_tot_bis_1=[T_tot_bis_1 T_1(:,i)];
    SD_tot_bis_1=[SD_tot_bis_1 ;SD_tot_1(i,:)];
    
    %i=i+sampling;
    i=i+1;
    j=j+1;
    else
        %i=i+sampling;
        i=i+1;
    end
end




    


E_FFSP_1=E_tot_bis_1;
SD_FFSP_1=SD_tot_bis_1;

%% 

f=figure;
f.Units='points';
f.OuterPosition=[10 10 1000 450];
freq=1;
subplot(1,2,1)
stairs(t_ob,Y(1,:),'m','LineWidth',2)
xlabel('t (s)')
ylabel('Y(t)')
xlim([0 t_ob(end)])
title('Observation Process')
set(gca,'FontSize',20)  


freq1=30;
subplot(1,2,2)
%stairs(t(t<=t_ob(end)),X_1_tot(1,t<=t_ob(end)),'b','LineWidth',2)
stairs(t,X_1_tot(1,:),'b','LineWidth',2)
hold on
errorbar(T_tot_bis([1 [1:freq:end]]),E_FFSP(1,[1 [1:freq:end]],1)',SD_FFSP(1,[1 [1:freq:end]])','r--','LineWidth',3)
hold on
errorbar(T_tot_bis_1([1 [1:freq:end]]),E_FFSP_1([1 [1:freq:end]],1)',SD_FFSP_1([1 [1:freq:end]],1)','--','Color',[0.4660 0.6740 0.1880],'LineWidth',3)
%errorbar(T_1([1 [1:freq1:end]]),E_tot_1([1 [1:freq1:end]],1)',SD_tot_1([1 [1:freq1:end]],1)','--','Color',[0.4660 0.6740 0.1880],'LineWidth',3)
% hold on
% errorbar(t_ob([1 [1:freq:end]]),E_pf([1 [1:freq:end]],2)',SD_pf([1 [1:freq:end]],1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold off
title('Hidden Process')
xlim([0 t_ob(end)])
%ylim([0 15])
xlabel('t (s)')
ylabel('Molecular Counts')
legend('Exact Trajectory', 'FFSP 1' ,'FFSP 2','Location','northwest')
set(gca,'FontSize',20)  


%%
%FilterDiff=abs(E-E_1);
ProbDiff=zeros(1,length(t_ob));
for i=1:length(t_ob)
    ProbDiff(i)=sum(abs(F(:,jump_times(i))-P(:,jump_times_1(i))));
end

Tab=[t_ob;eps;Err_jump_1;ProbDiff];
csvwrite('/Users/edambrosio/Desktop/PaperFigures/Toggle_Switch_Filter_1_Filter_2_bad_error',Tab);
