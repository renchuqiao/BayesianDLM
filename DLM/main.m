clear; 
kkk=10; % burn-in size

% EDA
[Ylog,Ynames,time]=xlsread('simulation_10genes.xlsx'); %Y is T x #_genomes i.e. T x 10
[T, q] = size(Ylog);
T = T-1;
q = q-1;
names=Ynames(1,2:q+1);



Y = (2.^Ylog)';
%setuptimeaxis  % to se
[q,T]=size(Y); 
k=10; delta=0.995; beta=0.975; r=0.99; n0=5; 
% k is size of data used to estimate prior 
% arp is lag 1

model=char('lag1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ma (Lt | Lt-1) for all genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 2; % N time lags in regression model
% PRED = zeros(q,T,N);
% RMSE = zeros(T,q,N);
% MAD = zeros(T,q,N);

% for arp = 1:N

arp = 1;
% lag = 1
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
%[a,b] = predictTVVAR(i,Yarp,Ynames,k,F,delta,beta,n0,r,arp);
w = i;
Y = Yarp;
T=size(Y,2);  
q=size(Y,1);  
my=mean(Y(:,1:k)')'; %take the first k days mean as empirical value for estimating the prior on intercept
vy=(std(Y(:,1:k)')').^2; %taking the standard dev of the first k days for estimating the prior on intercept 
p = size(F,1);  
h0=n0+q-1; %% What's h????
D0=h0*diag(vy);
z = zeros(p,q); zq=zeros(q,1); M0=z; 
        if (p>1)
              M0(2:1+q,:)=eye(q);     % sets prior mean to be zero everywhere but r on lag-1 of same series
              M0(1,:)=(1-r)*my; % and intercept accordingly based on 
        end
Mt = M0; C0=eye(p);  Ct=C0;  
n = n0; h=h0; D = D0;  St=D/h;
sMt=zeros(p,q,T);  sCt=zeros(p,p,T);  sdCt=zeros(p,q,T); sSt=zeros(q,q,T);  snt=zeros(1,T);  
sft=zeros(q,T); smft=sft; % FF 1-step forecasts, and BS fitted values
zq=zeros(q,1);sEt=zeros(q,T);sQt=zeros(q,q,T);At=zeros(q+arp,1);K=zeros(q,q);sloglik=zeros(1,T);

% Analysis -- 
% forward filtering: 
        for t = k+1:T
            ft = Mt'*F(:,t);            sft(:,t)=ft; 
            et = Y(:,t) - ft;
            Rt = Ct/delta; 
            h  = beta*h;  n=h-q+1;  D = beta*D;  snt(t)=n;     
            qvt = 1 + F(:,t)'*Rt*F(:,t); %qvt is a scaler
            sEt(:,t) = et./sqrt(qvt*diag(St)); %%%diag(St) is a 10x1 vector that returns the diagonal entries of matrix St
                                         sQt(:,:,t) = St*qvt; 
                                         sloglik(t) = ltpdf(et(w),0,qvt,n,D(w,w));
                                        
            At = Rt*F(:,t)/qvt;
            h=h+1; n=h-q+1; D = D+et*et'/qvt;  St=D/h; St=(St+St')/2;  
            Mt = Mt + At*et'; Ct = Rt - At*At'*qvt;   sCt(:,:,t)=Ct;
            sSt(:,:,t)=St; sMt(:,:,t) = Mt; sdCt(:,:,t) = sqrt(diag(Ct)*diag(St)'); 
        end
      
   %%Now Calculate log likelihood  
a=sloglik;
pred = zeros(1,T);
pred = sft(w,1:T);

% reverse smoothing  -
K=inv(sSt(:,:,T)); n=snt(T); Mt = sMt(:,:,T); Ct = sCt(:,:,T); smft(:,T) = Mt'*F(:,T);
for t=(T-1):-1:(k+1)
    n = (1-beta)*snt(t)+beta*snt(t+1);          snt(t)=n; 
    K=(1-beta)*inv(sSt(:,:,t))+beta*K;          St = inv(K); sSt(:,:,t)=St;  
    Mt = (1-delta)*sMt(:,:,t) +delta*Mt;        sMt(:,:,t) = Mt;   smft(:,t) = Mt'*F(:,t);
    Ct = (1-delta^2)*sCt(:,:,t) + delta^2*Ct;   sCt(:,:,t) = Ct; sdCt(:,:,t) = sqrt(diag(Ct)*diag(St)');
end            


b = pred;

l(i,:) = a;
mp(i,:) = b;
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% Explore some aspects of the retrospectivelly estimates of model quantities: 


% -1.  Some data, 1-step ahead forecasts and intervals  
%
k=10;     i=k+1:T; % plot only after first k 
for j=1:q
     ci = sqrt(squeeze(sQt(j,j,:))).*tinv(0.95,snt');  % credible interval width
     clf;
     figure(1); clf; 
     plot(i,Y(j,i),'rh',i,sft(j,i),'k--',i,sft(j,i)'+ci(i),'b:',i,sft(j,i)'-ci(i),'b:'); 
     legend('Data','Forecasts','Credible bands'); legend boxoff
     xlim([0 T+1]); %eval(xa)
     ylabel('Gene Expression Level');  title([names(j)]) 
     set(gca,'FontSize',20);
     figure(2); clf
     plot(i,Y(j,i)-sft(j,i),'k--',i,ci(i),'b:',i,-ci(i),'b:'); 
     legend('1-step errors','Credible bands'); legend boxoff
     xlim([0 T+1]); %eval(xa)
     ylabel('Prediction error for Gene Expression Level');  title([names(j)]) 
     set(gca,'FontSize',20);
     pause
end



% 2.  Some correlations over time 
%
figure(2); clf; 
for i = 1:q
    for j = (i+1):q
        J=[i,j]; 
        plot(1:T,squeeze(sSt(J(1),J(2),1:T))./sqrt(squeeze(sSt(J(1),J(1),1:T)).*squeeze(sSt(J(2),J(2),1:T))),'k')
        s = strcat('correlation between ', names(J(1)), ' and', names(J(2)))
        title(s)
        ylim([-1 1]); line([0 T+1],[0 0],'color','b','linestyle',':'); %eval(xa)
        set(gca,'FontSize',20);
        pause;
    end
end






