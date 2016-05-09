clear; 
kkk=10; % burn-in size


% EDA
[Ylog,Ynames,time]=xlsread('35toy.xlsx'); %Y is T x #_genomes i.e. T x 10
% original dataset genomicsLOG has gene expression level convted with log2
Ynom = char(Ynames);
[m_Ynom, n_Ynom] = size(Ynom);
Ynom = Ynom(49:48:m_Ynom,1:5);
% EDA plot on log2 scale
% plot_var = [Ylog(:,1),Ylog(:,2),Ylog(:,3),Ylog(:,4),Ylog(:,5),Ylog(:,6),Ylog(:,7),Ylog(:,8),Ylog(:,9),Ylog(:,10)];
% clf;
% plot(plot_var);
% legend(Ynom(1,:), Ynom(2,:), Ynom(3,:), Ynom(4,:), Ynom(5,:), Ynom(6,:), Ynom(7,:), Ynom(8,:), Ynom(9,:),Ynom(10,:))
% xlabel('Time points: 1 unit for 30 minutes, 24 hours in total');
% ylabel('log2 (gene expression level)');
% set(gca,'FontSize',20);

Y = (2.^Ylog)';
% convert data to original scale
 
%setuptimeaxis  % to se
[q,T]=size(Y); 
k=10; delta=0.995; beta=0.975; r=0.99; n0=5; 
% k is size of data used to estimate prior 
% arp is lag 1



model=char('lag1','lag2','lag3','lag4','lag5','lag6','lag7','lag8','lag9', ...
'lag10','lag11','lag12','lag13','lag14','lag15','lag16','lag17','lag18', ...
'lag20');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ma (Lt | Lt-1) for all genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 2; % N time lags in regression model
PRED = zeros(q,T,N);
RMSE = zeros(T,q,N);
MAD = zeros(T,q,N);

for arp = 1:N

F=zeros(q*arp,T-arp);
for j=1:arp, 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j));  
end
Yarp=Y; %Yarp=4xT
[q,T]=size(Yarp); 
Yarp(:,1:arp)=[];
F=[ones(1,T-arp);F];

l=zeros(q,T-arp); mp=zeros(q,T-arp);

for i=1:q
[a,b] = predictTVVAR(i,Yarp,Ynames,k,F,delta,beta,n0,r,arp);  
l(i,:)=a;
mp(i,:)=b;
end

meanpred=[zeros(q,arp) mp];    

pred_pointest = zeros(q,T);
rmse = zeros(T,q);
mad = zeros(T,q);

for w = 1:q
    pred_pointest(w,:) = meanpred(w,:);
    [rmse(:,w),mad(:,w)] = MSEMAD(w,pred_pointest(w,:),Y,kkk,T);
end

PRED(:,:,arp) = pred_pointest; 
RMSE(:,:,arp) = rmse;
MAD(:,:,arp) = mad;

end


% BestlagRMSE = zeros(q,2);
% BestlagMAD = zeros(q,2);
% 
% for w = 1:q
%     lastRMSE = [RMSE(T,w,1),RMSE(T,w,2),RMSE(T,w,3),RMSE(T,w,4),RMSE(T,w,5), ...
%     RMSE(T,w,6),RMSE(T,w,7),RMSE(T,w,8),RMSE(T,w,9),RMSE(T,w,10),RMSE(T,w,11), ...
%     RMSE(T,w,12),RMSE(T,w,13),RMSE(T,w,14),RMSE(T,w,15),RMSE(T,w,16),RMSE(T,w,17), ...
%     RMSE(T,w,18),RMSE(T,w,19),RMSE(T,w,20)];
% 
%     lastMAD = [MAD(T,w,1),MAD(T,w,2),MAD(T,w,3),MAD(T,w,4),MAD(T,w,5), ...
%     MAD(T,w,6),MAD(T,w,7),MAD(T,w,8),MAD(T,w,9),MAD(T,w,10),MAD(T,w,11), ...
%     MAD(T,w,12),MAD(T,w,13),MAD(T,w,14),MAD(T,w,15),MAD(T,w,16),MAD(T,w,17), ...
%     MAD(T,w,18),MAD(T,w,19),MAD(T,w,20)];
%     
%     [minRMSE,lagRMSE] = min(lastRMSE);
%     [minMAD,lagMAD] = min(lastMAD);
%     
%     BestlagRMSE(w,1) = lagRMSE;
%     BestlagRMSE(w,2) = minRMSE;
%     
%     BestlagMAD(w,1) = lagMAD;
%     BestlagMAD(w,2) = minMAD;
%     
%     plotRMSE = [RMSE(:,w,1),RMSE(:,w,2),RMSE(:,w,3),RMSE(:,w,4),RMSE(:,w,5), ...
%     RMSE(:,w,6),RMSE(:,w,7),RMSE(:,w,8),RMSE(:,w,9),RMSE(:,w,10),RMSE(:,w,11), ...
%     RMSE(:,w,12),RMSE(:,w,13),RMSE(:,w,14),RMSE(:,w,15),RMSE(:,w,16),RMSE(:,w,17), ...
%     RMSE(:,w,18),RMSE(:,w,19),RMSE(:,w,20)];
% 
%     plot(plotRMSE);
%     legend(model);
%     xlabel('Time points: 1 unit for 30 minutes, 24 hours in total');
%     ylabel(['Average Root Mean Squared Error']);
%     set(gca,'FontSize',20);
%     pause;
%     
%     
%     plotMAD = [MAD(:,w,1),MAD(:,w,2),MAD(:,w,3),MAD(:,w,4),MAD(:,w,5), ...
%     MAD(:,w,6),MAD(:,w,7),MAD(:,w,8),MAD(:,w,9),MAD(:,w,10),MAD(:,w,11), ...
%     MAD(:,w,12),MAD(:,w,13),MAD(:,w,14),MAD(:,w,15),MAD(:,w,16),MAD(:,w,17), ...
%     MAD(:,w,18),MAD(:,w,19),MAD(:,w,20)];
%     
%     
%     plot(plotMAD);
%     legend(model);
%     xlabel('Time points: 1 unit for 30 minutes, 24 hours in total');
%     ylabel(['Average Mean Absolute Deviation']);
%     set(gca,'FontSize',20);
%     pause;
%    
% end








