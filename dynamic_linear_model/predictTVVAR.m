function [l, pred] = predictTVVAR(w,Y,Ynames,k,F,delta,beta,n0,r,arp)
%%%%%%%%%%%%
%%%%before function, run [Y,Ynames,time]=xlsread('etfdata.xlsx');
%%%%%%% before function, specify F
%%%1) filename: data file name
%%%2) k = # of initial observations served for prior estimation; no FFBS
%%%3) F stores predictors
%%%4) delta=0.995; discount factor 
% discount model (level + TVVAR) 
% delta is the discount factor to estimate Vt (the evolutional variance
% matrix) where Wt = Pt / delta
% but what is Pt???? here in hte code, Rt = Ct / delta
% but what's Rt used for????
%%%5) beta =0.975;    % discount fator for volatility
% beta is used to est. the observational varaince matrix Vt
%Vt = (beta * Vt-1) / gamma_t where gamma_t 
%~ Be (beta * n_t-1 /2,(1-beta)*n_t-1 / 2)
%%%6) n0=5; (prior) n is the degrees of freedom in the final multivariate T distribution i.e. the n in the gamma distribution (gamma normal conjugacy)
%%%7) r=0.99; % priors %% r is the correlation between each two consecutive t
%%%8) arp: the lag number

 
T=size(Y,2);  
q=size(Y,1);  

% choose and fit a simple model: a simple local mean plus TV-VAR(1) model with volatility 
% first, transform to logs with a base of 2 ...
 
my=mean(Y(:,1:k)')'; %take the first k days mean as empirical value for estimating the prior on intercept
vy=(std(Y(:,1:k)')').^2; %taking the standard dev of the first k days for estimating the prior on intercept 

% lag is 1---- set up TV-VAR(arp) with locally constant level 

p = size(F,1);  % p is the dimension p is that of the number of predictors in F


%%%priors 
h0=n0+q-1; %% What's h????
D0=h0*diag(vy); %%create diagonal matrix q*q (10x10). Then fill in the diagon entries with the 10 entries in vy
%%What is the prior qxq (10x10) matrix D0 for?? Why set prior to be the
%%empirical variance?? WHAT IS D0???? 
z = zeros(p,q); zq=zeros(q,1); M0=z; 
        if (p>1)
              M0(2:1+q,:)=eye(q);     % sets prior mean to be zero everywhere but r on lag-1 of same series
              M0(1,:)=(1-r)*my; % and intercept accordingly based on 
        end
        % M_0 is prior on coefficients of all 
       
Mt = M0; C0=eye(p);  Ct=C0;        % initial Theta prior 
n = n0; h=h0; D = D0;  St=D/h;     % initial Sigma prior since D=h*diagon(vy), so St = diagon(vy)
 
% arrays to save model outputs 
sMt=zeros(p,q,T);  sCt=zeros(p,p,T);  sdCt=zeros(p,q,T); sSt=zeros(q,q,T);  snt=zeros(1,T);  
sft=zeros(q,T); smft=sft; % FF 1-step forecasts, and BS fitted values
zq=zeros(q,1);sEt=zeros(q,T);sQt=zeros(q,q,T);At=zeros(q+arp,1);K=zeros(q,q);sloglik=zeros(1,T);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
l=sloglik;
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

end