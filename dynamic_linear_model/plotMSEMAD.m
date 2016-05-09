function [plotMSE, plotMAD] = plotMSEMAD(mseSnew,mseSold,madSnew,madSold,new,old,w,Ynom,model,T,xa,tticks_3,tdates,kkk)

% Ynom(w,:) is SPY
% for example, new is 'M2', old is 'M1'

ratiomse=mseSnew./mseSold;
clf;
plotMSE=plot(1:T,ratiomse,'r');hline=refline([0 1]);
xlim([kkk-10 T]);
%axis([kkk-10,T,0.8,1.01])
eval(xa);set(gca,'FontSize',20);
title(['Ratio of MSE given by ',model(new,:),'over MSE given by ',model(old,:),'for daily closing prices of ',Ynom(w,:)],'FontSize', 24);
ylabel('Ratio of MSE','FontSize',36);

pause;


%%%%%%%%%%%%%%%% MAD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ratiomad=madSnew./madSold;
clf;
plotMAD=plot(1:T,ratiomad,'r');hline=refline([0 1]);
%axis([kkk-10,T,0.8,1.01])
xlim([kkk-10 T]);
eval(xa);set(gca,'FontSize',20);
title(['Ratio of MAD given by ',model(new,:),'over MAD given by ',model(old,:),'for daily closing prices of ',Ynom(w,:)],'FontSize', 24);
ylabel('Ratio of MAD','FontSize',36);

pause;

