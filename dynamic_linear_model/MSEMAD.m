function [mseS, madS] = MSEMAD(w,meanpred,Y_nonfx,kkk,T)

% predicted gene level actual scale:
pred=meanpred;
% realized gene level actual scale:
level=Y_nonfx(w,:);

% levels MSE and MAD
diff=level-pred;  %error term. Where diff dimension is qxT 
diff(1:kkk)=[];
diff=[zeros(1,kkk) diff];    % now diff is 1xT

sqeS=zeros(1,T);   

for t=1:T,
    e=diff(t);   
    sqeS(t)=(e^2);
end

mseS=zeros(1,T);
for t=1:T,
    mseS(t)=sqrt((sum(sqeS(1:t)))/t);  % w is the # of the stock, w=1 is GLD; w=2 is QQQ; w=3 is SPY; w=4 is USO    
end


%   Root-mean-square error http://en.wikipedia.org/wiki/Root-mean-square_deviation




% MAD
madS=zeros(1,T);
diffabs=abs(diff); % here diff is 1xT dimension 

for t=1:T,
    madS(t)=(sum(diffabs(1:t)))/t;   % w is the # of the stock, w=1 is GLD; w=2 is QQQ; w=3 is SPY; w=4 is USO
end




