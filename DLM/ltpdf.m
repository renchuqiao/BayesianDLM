function ltpdf = ltpdf(x,m,w,n,D)
% q-dim marginal t pdf evaluated at x:  
% x model:      x|Sigma ~ N( m, w.Sigma)   
% CONJ PRIOR:     Sigma ~ IW_p(n,D): n>0, d=n+q-1, d>q-1, D=prior sumofsqs
%                                               est Sigma=D/n = S
% output: pdf = p(x) 
%
q=size(x,1); C=chol(D)'; e=inv(C)*(x-m);  d=n+q-1; 
ltpdf = q*log(2)/2 - q*log(2*pi*w)/2 -sum(log(diag(C))) -(d+1)*log(1+(e'*e)/w)/2;
ltpdf = ltpdf + sum(gammaln((1+d-(0:q-1))/2) - gammaln((d-(0:q-1))/2)); 
 
    

