function [results] =mcmcwfbd8b(nit,br, n,T,Y,YL,nr, X1,ew,B1,BL,E1, Fw,Fw19L,M,MM,W5,W5L,W5g,fiv,inpop,wea1,wea2,wea3,wea4,lambda0, rho0,mu0,gamma0,lgamma,ugamma,beta0, delta0,phi0,eta0,kappa0,zeta0,sigmav0,sigmae0,sigmat0,sigmam0,X47,indiv10,indiv20,indiv30,Time10,Time30,rs,sce)
% =========================================================================
% This function conducts the MCMC estimation of a higher-order SAR model with spatial errors.
%==========================================================================

nomit=nit*br; %burn-in sample
   
%store all parameters
% 
% popweight=X47./(sum(X47));

%parameter in the SDPD equation
lambdas=zeros(nit,1);  gamma1s=zeros(nit,1); gamma2s=zeros(nit,1);
 rhods=zeros(nit,1);
rho1s=zeros(nit,1); rho2s=zeros(nit,1); rho3s=zeros(nit,1); rho4s=zeros(nit,1); mus=zeros(nit,1);  beta1s=zeros(nit,1); beta2s=zeros(nit,1); beta3s=zeros(nit,1); beta4s=zeros(nit,1); beta5s=zeros(nit,1);
beta6s=zeros(nit,1); sigmavs=zeros(nit,1); kappa1s=zeros(nit,1); zeta1s=zeros(nit,1);   betaps=zeros(nit,1); sigmacs=zeros(nit,1);
 
%parameter for the linear equation on vaccination
delta1s=zeros(nit,1);  delta2s=zeros(nit,1); delta3s=zeros(nit,1); delta4s=zeros(nit,1); delta5s=zeros(nit,1); delta6s=zeros(nit,1);
sigmaes=zeros(nit,1);  kappa2s=zeros(nit,1); indiv2s=zeros(n,1);  zeta2s=zeros(nit,1);   

%parameter for the linear equation on within state flow
phi1s=zeros(nit,1);  phi2s=zeros(nit,1); phi3s=zeros(nit,1); phi4s=zeros(nit,1); phi5s=zeros(nit,1); phi6s=zeros(nit,1);
kappa3s=zeros(nit,1); sigmats=zeros(nit,1);  phips=zeros(nit,1);  zeta3s=zeros(nit,1); sigmaps=zeros(nit,1);
  
  
%parameter for gravity equation
kappa4s=zeros(nit,1);  zeta4s=zeros(nit,1); eta1s=zeros(nit,1); eta2s=zeros(nit,1); eta3s=zeros(nit,1); eta4s=zeros(nit,1); eta5s=zeros(nit,1);
eta6s=zeros(nit,1); eta7s=zeros(nit,1); eta8s=zeros(nit,1); eta9s=zeros(nit,1); eta10s=zeros(nit,1);
sigmams=zeros(nit,1);

%store transform data
Yr=zeros(n,T); D0=cell(T,1); M0=cell(T,1);  D1=cell(T,1); M1=cell(T,1);
Yc=zeros(T,nit); % baseline counterfactual national new case

Wp=zeros(n,n);  stateoutflow=zeros(n,T,nit-nomit); stateinflow=zeros(n,T,nit-nomit);
outflowm=zeros(n,T); outflowu=zeros(n,T); outflowl=zeros(n,T);
inflowm=zeros(n,T); inflowu=zeros(n,T); inflowl=zeros(n,T);


%store national CF new case
Ymean=zeros(n,T);
Yub=zeros(n,T);  Ylb=zeros(n,T);

%store national CF cumulative case
Ycmean=zeros(n,1);
Ycub=zeros(n,1);  Yclb=zeros(n,1);

stateccase=kron(YL(:,1),ones(1,nit-nomit)); % initial CF cumulative case 

Vc=zeros(T,nit); %store CF posterior mean of national average vaccination.


Vmean=zeros(n,T);
Vub=zeros(n,T);  Vlb=zeros(n,T);


statecase=zeros(n,T,nit-nomit);
statevax=zeros(n,T,nit-nomit);

%=======================================
%prepare for the AM algorithm for lambda 
%=======================================
KKt=4;   st1=0.1/sqrt(KKt); ratio=0.05;
sum1=zeros(6,1);  sum2=zeros(6,6);
theta0=[lambda0;rho0;mu0];

KKg=2; sg1=0.1/sqrt(KKg); sumg1=zeros(2,1); sumg2=zeros(2,2);

%==========================
% value of prior parameters
%===========================
TT=200; %initial period of AM algorithm 
P=1000000;  % for beta1 and beta2   
a=0.001;b=0.001;  %for sigma
c=0.001; d=0.001;

sigmac0=1; betap0=1;
sigmap0=1; phip0=1;

omega0=cell(T,1); xi0=cell(T,1);
for t=1:T
   omega0{t}=zeros(n,1); 
    xi0{t}=zeros(n,1);
end

kappa10=kappa0(1); kappa20=kappa0(2); kappa30=kappa0(3); kappa40=kappa0(4);
zeta10=zeta0(1); zeta20=zeta0(2); zeta30=zeta0(3); zeta40=zeta0(4);



% Xd47=kron(X47,ones(T,1));
% Xo47=kron(ones(T,1),X47);


BBL=MM*BL;
inpopr=eye(n)-inpop;

for t=1:T
      D0{t}=diag(B1(:,t)<=gamma0(1));
      M0{t}=diag(B1(:,t)<=gamma0(2));
  YL1=D0{t}*inpop*YL(:,t); YL2=(eye(n)-D0{t})*inpop*YL(:,t);
  YL3=M0{t}*inpopr*YL(:,t);
  YL4=(eye(n)-M0{t})*inpopr*YL(:,t);
    
  W10=W5{t}; 
   

    SS=eye(n)-lambda0*W10;
 Yr(:,t)=SS*Y(:,t)-rho0(1)*YL1-rho0(2)*YL2-rho0(3)*YL3-rho0(4)*YL4-mu0*W5L{t}*YL(:,t); 
    
end

for i=1: nit
    
     if mod(i,10)==0
        fprintf('%d\n',i)
         fprintf('lambda0:   %.3f %.3f %.3f\n', lambda0);
          fprintf('rho10:   %.3f %.3f %.3f\n', rho0(1));
          fprintf('rho20:   %.3f %.3f %.3f\n', rho0(2));
          fprintf('rho30:   %.3f %.3f %.3f\n', rho0(3));
          fprintf('rho40:   %.3f %.3f %.3f\n', rho0(4));
          fprintf('gamma10: %.3f %.3f %.3f\n', gamma0(1));
             fprintf('gamma20: %.3f %.3f %.3f\n', gamma0(2));
%              
%                fprintf('beta10: %.3f %.3f %.3f\n', beta0(1));
          

           fprintf('delta10: %.3f %.3f %.3f\n', delta0(1));
              fprintf('delta20: %.3f %.3f %.3f\n', delta0(2));
            fprintf('phi10: %.3f %.3f %.3f\n', phi0(1));
             fprintf('phi20: %.3f %.3f %.3f\n', phi0(2));
%                 
%              fprintf('eta10: %.3f %.3f %.3f\n', eta0(1));
             fprintf('eta20: %.3f %.3f %.3f\n', eta0(2));
              fprintf('eta30: %.3f %.3f %.3f\n', eta0(3));
%               fprintf('sigmam0: %.3f %.3f %.3f\n', sigmam0);


     end
     
     
     % sample destination latent omega
N=n-1;
sigma0i=[1/sigmav0, 0, 0,zeros(1,N) ; 0, 1/sigmae0,0, zeros(1,N) ; 0, 0, 1/sigmat0,zeros(1,N); zeros(N,1),zeros(N,1),zeros(N,1),eye(N)/sigmam0]; kappar0=[kappa10; kappa20;kappa30;ones(N,1)*kappa40];
 %kappab=kron(kappar0,eye(n)); 
  kappab=kron(eye(n),kappar0);
 %Sigma0i=kron(sigma0i,eye(n));
 Sigma0i=kron(eye(n),sigma0i);
 Sigmaw=(eye(n)+kappab'*Sigma0i*kappab)\eye(n);
 
for t=1:T
  YY=zeros((N+3)*n,1); fivt=fiv{t}; w5gt=W5g{t}; xit=xi0{t};
  wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
 
  
   for ii=1:n
    yy4=[]; 
     zz1=[Fw(ii,t), X1(ii,:,t),nr(ii,:,t)];
   
     zz2=[BBL(ii,t),E1(ii,t),ew(ii,:,t)];
     
     zz3=[B1(ii,t),Fw19L(ii,t), X1(ii,:,t)];
     
     yy1=Yr(ii,t)-zz1*beta0-Time10(t)-indiv10(ii)-xit(ii)*zeta10;
     yy2=B1(ii,t)-zz2*delta0-indiv20(ii)-xit(ii)*zeta20;
     yy3=Fw(ii,t)-zz3*phi0-indiv30(ii)-Time30(t)-xit(ii)*zeta30;

     for jj=1:n
        if (ii~=jj)
          r=[fivt(ii,jj),B1(ii,t),B1(jj,t),wea1t(ii,jj),wea2t(ii,jj),wea3t(ii,jj),wea4t(ii,jj),X47(ii),X47(jj),M(ii,jj)]; 
        gg=w5gt(ii,jj)-r*eta0-xit(jj)*zeta40;
         yy4=[yy4;gg];
        end
     end
     
     




     
     yy=[yy1;yy2;yy3;yy4];
     YY((ii-1)*51+1:ii*51)=yy;

   end

   Tw=Sigmaw*kappab'*(Sigma0i)*YY;
   omegaft=mvnrnd(Tw,Sigmaw);
   omegafn=omegaft';

   omega0{t}=omegafn;
 
end




% 
% for t=1:T
%   %  omega0{t}=omegas(:,t);
%     if (i>nomit)
%     omegasum(:,t)=omegasum(:,t)+omega0{t};
%     
%     end
% end


     % sample origin latent omega
N=n-1;
sigma0i=[1/sigmav0, 0, 0,zeros(1,N) ; 0, 1/sigmae0,0, zeros(1,N) ; 0, 0, 1/sigmat0,zeros(1,N); zeros(N,1),zeros(N,1),zeros(N,1),eye(N)/sigmam0]; zetar0=[zeta10; zeta20;zeta30;ones(N,1)*zeta40];
 %kappab=kron(kappar0,eye(n)); 
  zetab=kron(eye(n),zetar0);
 %Sigma0i=kron(sigma0i,eye(n));
 Sigma0i=kron(eye(n),sigma0i);
 Sigmax=(eye(n)+zetab'*Sigma0i*zetab)\eye(n);
 
for t=1:T
  YY=zeros((N+3)*n,1); fivt=fiv{t}; w5gt=W5g{t}; omegat=omega0{t};
  wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
  
   for ii=1:n
    yy4=[]; 
     zz1=[Fw(ii,t), X1(ii,:,t),nr(ii,:,t)];
   
     zz2=[BBL(ii,t),E1(ii,t),ew(ii,:,t)];
     
     zz3=[B1(ii,t),Fw19L(ii,t),X1(ii,:,t)];
     
     yy1=Yr(ii,t)-zz1*beta0-Time10(t)-indiv10(ii)-omegat(ii)*kappa10;
     yy2=B1(ii,t)-zz2*delta0-indiv20(ii)-omegat(ii)*kappa20;
     yy3=Fw(ii,t)-zz3*phi0-Time30(t)-indiv30(ii)-omegat(ii)*kappa30;

     for jj=1:n
        if (ii~=jj)
          r=[fivt(jj,ii),B1(jj,t),B1(ii,t),wea1t(jj,ii),wea2t(jj,ii),wea3t(jj,ii),wea4t(jj,ii),X47(jj),X47(ii),M(jj,ii)]; 
        gg=w5gt(jj,ii)-r*eta0-omegat(jj)*kappa40;
         yy4=[yy4;gg];
        end
     end
     
     




     
     yy=[yy1;yy2;yy3;yy4];
     YY((ii-1)*51+1:ii*51)=yy;

   end

   Tx=Sigmax*zetab'*(Sigma0i)*YY;
   xift=mvnrnd(Tx,Sigmax);
   xifn=xift';

   xi0{t}=xifn;
 
end





% for t=1:T
%   %  omega0{t}=omegas(:,t);
%     if (i>nomit)
%     xisum(:,t)=xisum(:,t)+xi0{t};
%     
%     end
% end
   
     
     
     
     
    
   if (i<=2*TT)
      accept=0;
      theta1=mvnrnd(theta0,st1^2*eye(6));
      theta1=theta1';
   lambda1=theta1(1);
      rho1=theta1(2:5);
      mu1=theta1(6);
  
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
     
     
      
      while (accept==0) %reject bounds on lambda1 and rho1
          if (lambdam*rs<1 ) && (lambdam*rs+rhom+mum*rs<1)
              accept=1;
          else
           theta1=mvnrnd(theta0,st1^2*eye(6));
      theta1=theta1';
   lambda1=theta1(1);
      rho1=theta1(2:5);
      mu1=theta1(6);

     
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
     
          end
      end
  end
  
  if (i>2*TT)
      accept=0;
      if (i<=nomit)
      vvarr=(1-ratio)^2*2.38^2*varr/KKt+ratio^2*st1^2*eye(6);
      end
      theta1=mvnrnd(theta0,vvarr);
      theta1=theta1';
 
  
   lambda1=theta1(1);
      rho1=theta1(2:5);
      mu1=theta1(6);
     
     
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
      
      while (accept==0) %reject bounds on lambda1
          if (lambdam*rs<1 ) && (lambdam*rs+rhom+mum*rs<1)
              accept=1;
          else
      theta1=mvnrnd(theta0,vvarr);
    theta1=theta1';
 
    lambda1=theta1(1);
      rho1=theta1(2:5);
      mu1=theta1(6);
     
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
      mum=max(abs(mu1));
     
          end
      end
  end

   
 
 frpp=1;
 

  
  for t=1: T
    
      D0{t}=diag(B1(:,t)<=gamma0(1));
      M0{t}=diag(B1(:,t)<=gamma0(2));
  YL10=D0{t}*inpop*YL(:,t); YL20=(eye(n)-D0{t})*inpop*YL(:,t);
  YL30=M0{t}*inpopr*YL(:,t);
  YL40=(eye(n)-M0{t})*inpopr*YL(:,t);  
      
      



 
       W10=W5{t}; 
   
   
      
 frpp=likn(Y(:,t), YL(:,t), YL10,YL20,YL30,YL40, YL10,YL20, YL30, YL40, nr(:,:,t),X1(:,:,t), Fw(:,t), W10,W5L{t},omega0{t},xi0{t},indiv10, Time10(t),lambda0,rho0,mu0,lambda1,rho1,mu1,beta0,kappa10,zeta10,sigmav0)*frpp;
  end
    

    
   % Determine the transition probability
    Acceptr=min(1,frpp);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1) 
        lambda_1=lambda1;
        rho_1=rho1;
        mu_1=mu1;
 
       
     else
       lambda_1=lambda0;
        rho_1=rho0;
        mu_1=mu0;
  
       
         
     end
    
    % Store the value and continue
    rho0=rho_1;
  rho1s(i,1)=rho0(1);
  rho2s(i,1)=rho0(2);
  rhods(i,1)=rho0(1)-rho0(2);
  rho3s(i,1)=rho0(3);
  rho4s(i,1)=rho0(4);
       
     mu0=mu_1;
     mus(i,1)=mu0;
     
  
     lambda0=lambda_1;
     lambdas(i,1)=lambda0;
    

      

      
  %calculate the empirical covariance
    theta0=[lambda0;rho0; mu0];
    
    if (i<=nomit)
        sum1=sum1+theta0;
        sum2=sum2+theta0*theta0';
    end
    
    if (i>1)&&(i<=nomit)
        mean1=sum1/i;
        varr=sum2/i-mean1*mean1';
    end
 
    
    %sample gamma
    
    if (i<=2*TT)
        accept=0;
        gamma1=mvnrnd(gamma0,sg1^2*eye(2));
        
        while (accept==0) %reject bounds on lambda1 and rho1
            if (gamma1(1)<ugamma ) && (gamma1(1)>lgamma)&& (gamma1(2)<ugamma ) && (gamma1(2)>lgamma)
                accept=1;
            else
                gamma1=mvnrnd(gamma0,sg1^2*eye(2));
                
            end
        end
    end
    
    if (i>2*TT)
        accept=0;
        if (i<=nomit)
            vvarrg=(1-ratio)^2*2.38^2*varrg/KKg+ratio^2*sg1^2;
        end
        gamma1=mvnrnd(gamma0,vvarrg);
        
        
        while (accept==0) %reject bounds on lambda1
            if (gamma1(1)<ugamma ) && (gamma1(1)>lgamma)&& (gamma1(2)<ugamma ) && (gamma1(2)>lgamma)
                accept=1;
            else
                gamma1=mvnrnd(gamma0,vvarrg);
                
            end
        end
    end
    
    
    frpp=1;
    
    for t=1: T
        D0{t}=diag(B1(:,t)<=gamma0(1));
        D1{t}=diag(B1(:,t)<=gamma1(1));
        
        M0{t}=diag(B1(:,t)<=gamma0(2));
        M1{t}=diag(B1(:,t)<=gamma1(2));
        
        W10=W5{t}; 
        
        YL10=D0{t}*inpop*YL(:,t); YL20=(eye(n)-D0{t})*inpop*YL(:,t);
        YL30=M0{t}*inpopr*YL(:,t); YL40=(eye(n)-M0{t})*inpopr*YL(:,t);
     
        
         YL11=D1{t}*inpop*YL(:,t); YL21=(eye(n)-D1{t})*inpop*YL(:,t);
        YL31=M1{t}*inpopr*YL(:,t); YL41=(eye(n)-M1{t})*inpopr*YL(:,t);
        
   
        
        frpp=likn(Y(:,t), YL(:,t), YL10,YL20,YL30,YL40,YL11,YL21,YL31,YL41, nr(:,:,t), X1(:,:,t),Fw(:,t), W10,W5L{t},omega0{t},xi0{t}, indiv10,  Time10(t),lambda0,rho0,mu0,lambda0,rho0,mu0,beta0,kappa10,zeta10,sigmav0)*frpp;
    end
    
    
    
    % Determine the transition probability
    Acceptr=min(1,frpp);
    
    % Draw from uniform(0,1)
    u1=rand(1,1);
    
    % Transition to candidate delta11 with probability acceptr
    if (Acceptr>u1)
        gamma_1=gamma1;
    else
        gamma_1=gamma0;
    end
    
    % Store the value and continue
    gamma0=gamma_1;
    gamma1s(i,1)=gamma0(1);
    gamma2s(i,1)=gamma0(2);
    
    
    
    
    
    %calculate the empirical covariance
    
    
    if (i<=nomit)
        sumg1=sumg1+gamma0;
        sumg2=sumg2+gamma0*gamma0';
    end
    
    if (i>1)&&(i<=nomit)
        mean1=sumg1/i;
        varrg=sumg2/i-mean1*mean1';
    end
    
    
    
    
      
for t=1:T
      D0{t}=diag(B1(:,t)<=gamma0(1));
      M0{t}=diag(B1(:,t)<=gamma0(2));
      
  YL1=D0{t}*inpop*YL(:,t); YL2=(eye(n)-D0{t})*inpop*YL(:,t); 
   YL3=M0{t}*inpopr*YL(:,t); YL4=(eye(n)-M0{t})*inpopr*YL(:,t); 
     W10=W5{t}; 
   

    SS=eye(n)-lambda0(1)*W10;
 Yr(:,t)=SS*Y(:,t)-rho0(1)*YL1-rho0(2)*YL2-rho0(3)*YL3-rho0(4)*YL4-mu0*W5L{t}*YL(:,t); 
    
end
    

   


  %sample beta0=[beta10; beta20;beta30;beta40;betar0]
 sumb=zeros(6,1); sumbv=zeros(6,6);
 
for t=1:T
 ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
   omegat=omega0{t}; xit=xi0{t};
   sumb=sumb+ZZ'*(Yr(:,t)-indiv10-ones(n,1)*Time10(t)-omegat*kappa10-xit*zeta10)/sigmav0;
 sumbv=sumbv+ZZ'*ZZ/sigmav0;

end

sumv=(sumbv+eye(6)/(P))\eye(6);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumb;

betat=mvnrnd(Tr,sumv); 


beta0=betat';


beta10=beta0(1);
beta20=beta0(2);
beta30=beta0(3);
beta40=beta0(4);
beta50=beta0(5);
beta60=beta0(6);

beta1s(i,1)=beta10;
beta2s(i,1)=beta20;
beta3s(i,1)=beta30;
beta4s(i,1)=beta40;
beta5s(i,1)=beta50;
beta6s(i,1)=beta60;


%Sample kappa10
sumk=0; sumkv=0;
for t=1:T
ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
  omegat=omega0{t};  xit=xi0{t};
sumk=sumk+omegat'*(Yr(:,t)-indiv10-ones(n,1)*Time10(t)-ZZ*beta0-xit*zeta10)/sigmav0;
 sumkv=sumkv+omegat'*omegat/sigmav0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;

kappa10=mvnrnd(Tr,sumv); 

kappa1s(i,1)=kappa10;

%Sample zeta10
sumk=0; sumkv=0;
for t=1:T
ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
  omegat=omega0{t};  xit=xi0{t};
sumk=sumk+xit'*(Yr(:,t)-indiv10-ones(n,1)*Time10(t)-ZZ*beta0-omegat*kappa10)/sigmav0;
 sumkv=sumkv+xit'*xit/sigmav0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;

zeta10=mvnrnd(Tr,sumv); 

zeta1s(i,1)=zeta10;




%Sample sigmav
 sums=0;
 for t=1:T
 ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
    omegat=omega0{t}; xit=xi0{t};
  ss=Yr(:,t)-ZZ*beta0-indiv10-ones(n,1)*Time10(t)-omegat*kappa10-xit*zeta10;
  sums=sums+ss'*ss;   
 end
 
 
 ap=a+n*T/2;
 bp=b+0.5*sums;
sigmav0= 1./gamrnd(ap,1/bp);  
sigmavs(i,1)=sigmav0;   


 %Gibbs sampling for time-fixed effects
    
    for t=2:T
      ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
          omegat=omega0{t}; xit=xi0{t};
        Yf=Yr(:,t)-ZZ*beta0-indiv10-omegat*kappa10-xit*zeta10;
        
        Sigmaf=(0.0001+n/(sigmav0))^(-1);
       
        Tf=Sigmaf*ones(n,1)'*Yf/(sigmav0);
        
        time0=mvnrnd(Tf,Sigmaf);
        
        Time10(t,1)=time0;
       
    end
    




    
 %sample individual fixed effects C1
sumr=zeros(n,n); summ=zeros(n,1);
for q2=1:T
sumr=sumr+eye(n)/sigmav0;
ZZ=[Fw(:,q2), X1(:,:,q2),nr(:,:,q2)];
  omegat=omega0{q2}; xit=xi0{q2};
Tc=Yr(:,q2)-ZZ*beta0-ones(n,1)*Time10(q2)-omegat*kappa10-xit*zeta10;
summ=summ+Tc/sigmav0;
end
Sigmai=(eye(n)/sigmac0+sumr)\eye(n);
%Sigmai=mean(cat(3,Sigmai,Sigmai'),3);   
Ti=Sigmai*(summ+X47*betap0/sigmac0);
indiv11=mvnrnd(Ti,Sigmai); 
indiv10=indiv11';




%sample betap0
Sigmai=(0.0001+X47'*X47/sigmac0)^(-1);
Ti=Sigmai*(X47'*indiv10/sigmac0);
betap0=mvnrnd(Ti,Sigmai);
betaps(i,1)=betap0;


%Sample sigmac

 
 
 cp=c+n/2;
 dp=d+0.5*(indiv10-X47*betap0)'*(indiv10-X47*betap0);
sigmac0= 1./gamrnd(cp,1/dp);  
sigmacs(i,1)=sigmac0;  




%*****************************
% Sample parameters in the linear panel on vaccination
%******************************
    sume=zeros(6,1); sumev=zeros(6,6);
 
for t=1:T
 ZZ=[BBL(:,t),E1(:,t),ew(:,:,t)];
   omegat=omega0{t}; xit=xi0{t};
 sume=sume+ZZ'*(B1(:,t)-indiv20-omegat*kappa20-xit*zeta20)/sigmae0;
 sumev=sumev+ZZ'*ZZ/sigmae0;

end

sumv=(sumev+eye(6)/(P))\eye(6);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sume;

delta0=mvnrnd(Tr,sumv); 
delta0=delta0';

delta1s(i,1)=delta0(1);
delta2s(i,1)=delta0(2);
delta3s(i,1)=delta0(3);
delta4s(i,1)=delta0(4);
delta5s(i,1)=delta0(5);
delta6s(i,1)=delta0(6);

%Sample sigmae
 sums=0;
 for t=1:T
    ZZ=[BBL(:,t),E1(:,t),ew(:,:,t)];
     omegat=omega0{t}; xit=xi0{t}; 
 ss=B1(:,t)-ZZ*delta0-indiv20-omegat*kappa20-xit*zeta20;
  sums=sums+ss'*ss;   
 end
 
 
 ap=a+n*T/2;
 bp=b+0.5*sums;
sigmae0= 1./gamrnd(ap,1/bp);  
sigmaes(i,1)=sigmae0;   

 

 % sample fixed effect
 sumr=zeros(n,n); summ=zeros(n,1);
for q2=1:T
sumr=sumr+eye(n)/sigmae0;
 ZZ=[BBL(:,q2),E1(:,q2),ew(:,:,q2)];
 omegat=omega0{q2}; xit=xi0{q2};
Tc=B1(:,q2)-ZZ*delta0-omegat*kappa20-xit*zeta20;
summ=summ+Tc/sigmae0;
end
Sigmai=(0.0001*eye(n)+sumr)\eye(n);
%Sigmai=mean(cat(3,Sigmai,Sigmai'),3);   
Ti=Sigmai*summ;
indiv21=mvnrnd(Ti,Sigmai); 
indiv20=indiv21';

if (i>nomit)
indiv2s=indiv2s+indiv20;
end


 
% sample kappa20
sumk=0; sumkv=0;
for t=1:T
     ZZ=[BBL(:,t),E1(:,t),ew(:,:,t)];
  omegat=omega0{t}; xit=xi0{t};
  sumk=sumk+omegat'*(B1(:,t)-indiv20-ZZ*delta0-xit*zeta20)/sigmae0;
 sumkv=sumkv+omegat'*omegat/sigmae0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumk;


 kappa20=mvnrnd(Tr,sumv); 

kappa2s(i,1)=kappa20;

% sample zeta20
sumk=0; sumkv=0;
for t=1:T
     ZZ=[BBL(:,t),E1(:,t),ew(:,:,t)];
  omegat=omega0{t}; xit=xi0{t};
  sumk=sumk+xit'*(B1(:,t)-indiv20-ZZ*delta0-omegat*kappa20)/sigmae0;
 sumkv=sumkv+xit'*xit/sigmae0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumk;


 zeta20=mvnrnd(Tr,sumv); 

zeta2s(i,1)=zeta20;


  


%*****************************
% Sample parameters in the linear panel on within state flow
%******************************
    sumt=zeros(6,1); sumtv=zeros(6,6);
 
for t=1:T
 ZZ=[B1(:,t),Fw19L(:,t),X1(:,:,t)];
   omegat=omega0{t}; xit=xi0{t};
 sumt=sumt+ZZ'*(Fw(:,t)-indiv30-omegat*kappa30-xit*zeta30-ones(n,1)*Time30(t))/sigmat0;
 sumtv=sumtv+ZZ'*ZZ/sigmat0;

end

sumv=(sumtv+eye(6)/(P))\eye(6);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumt;

phi0=mvnrnd(Tr,sumv); 
phi0=phi0';

phi1s(i,1)=phi0(1);
phi2s(i,1)=phi0(2);
phi3s(i,1)=phi0(3);
phi4s(i,1)=phi0(4);
phi5s(i,1)=phi0(5);
phi6s(i,1)=phi0(6);

%Sample sigmat
 sums=0;
 for t=1:T
      ZZ=[B1(:,t),Fw19L(:,t),X1(:,:,t)];
     omegat=omega0{t}; xit=xi0{t};
 ss=Fw(:,t)-indiv30-ZZ*phi0-omegat*kappa30-xit*zeta30-ones(n,1)*Time30(t);
  sums=sums+ss'*ss;   
 end
 
 
 ap=a+n*T/2;
 bp=b+0.5*sums;
sigmat0= 1./gamrnd(ap,1/bp);  
sigmats(i,1)=sigmat0; 

 %Gibbs sampling for time-fixed effects
    
    for t=2:T
      ZZ=[B1(:,t),Fw19L(:,t),X1(:,:,t)];
          omegat=omega0{t};xit=xi0{t};
        Yf=Fw(:,t)-indiv30-ZZ*phi0-omegat*kappa30-xit*zeta30;
        
        Sigmaf=(0.0001+n/(sigmat0))^(-1);
       
        Tf=Sigmaf*ones(n,1)'*Yf/(sigmat0);
        
        time0=mvnrnd(Tf,Sigmaf);
        
        Time30(t,1)=time0;
       
    end
  






%  
% sample kappa30
sumk=0; sumkv=0;
for t=1:T
     ZZ=[B1(:,t),Fw19L(:,t),X1(:,:,t)];
  omegat=omega0{t}; xit=xi0{t};
sumk=sumk+omegat'*(Fw(:,t)-ZZ*phi0-indiv30-ones(n,1)*Time30(t)-xit*zeta30)/sigmat0;
 sumkv=sumkv+omegat'*omegat/sigmat0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumk;

% kappa20=normt_rnd(Tr,sumv,-9999999,0);
 kappa30=mvnrnd(Tr,sumv); 

kappa3s(i,1)=kappa30;


% sample zeta30
sumk=0; sumkv=0;
for t=1:T
     ZZ=[B1(:,t),Fw19L(:,t),X1(:,:,t)];
  omegat=omega0{t}; xit=xi0{t};
sumk=sumk+xit'*(Fw(:,t)-ZZ*phi0-indiv30-ones(n,1)*Time30(t)-omegat*kappa30)/sigmat0;
 sumkv=sumkv+xit'*xit/sigmat0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumk;

% kappa20=normt_rnd(Tr,sumv,-9999999,0);
 zeta30=mvnrnd(Tr,sumv); 

zeta3s(i,1)=zeta30;


% sample fixed effect
 sumr=zeros(n,n); summ=zeros(n,1);
for q2=1:T
sumr=sumr+eye(n)/sigmat0;
 ZZ=[B1(:,q2),Fw19L(:,q2),X1(:,:,q2)];
 omegat=omega0{q2}; xit=xi0{q2};
Tc=Fw(:,q2)-ZZ*phi0-omegat*kappa30-xit*zeta30-ones(n,1)*Time30(q2);
summ=summ+Tc/sigmat0;
end
Sigmai=(eye(n)/sigmap0+sumr)\eye(n);
%Sigmai=mean(cat(3,Sigmai,Sigmai'),3);   
Ti=Sigmai*(summ+X47*phip0/sigmap0);
indiv31=mvnrnd(Ti,Sigmai); 
indiv30=indiv31';



%sample phip0
Sigmai=(0.0001+X47'*X47/sigmap0)^(-1);
Ti=Sigmai*(X47'*indiv30/sigmap0);
phip0=mvnrnd(Ti,Sigmai);
phips(i,1)=phip0;


%Sample sigmap
cp=c+n/2;
 dp=d+0.5*(indiv30-X47*phip0)'*(indiv30-X47*phip0);
sigmap0= 1./gamrnd(cp,1/dp);  
sigmaps(i,1)=sigmap0;  






%**************************
% gravity equation
%*************************

sumeta=zeros(10,1); sumetav=zeros(10,10);
for t=1:T

    fivt=fiv{t}; omegat=omega0{t};  w5gt=W5g{t}; xit=xi0{t};
    wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
  
    for ii=1:n
        for jj=1:n
         
         
            
        if (ii~=jj)
         r=[fivt(ii,jj),B1(ii,t),B1(jj,t),wea1t(ii,jj),wea2t(ii,jj),wea3t(ii,jj),wea4t(ii,jj),X47(ii),X47(jj),M(ii,jj)]; 
        sumeta=sumeta+r'*(w5gt(ii,jj)-omegat(ii)*kappa40-xit(jj)*zeta40)/sigmam0;
        sumetav=sumetav+r'*r/sigmam0;
        
        end
        
        

        end
    end
end

sumv=(sumetav+eye(10)/(P))\eye(10);

Tr=sumv*sumeta;

eta0=mvnrnd(Tr,sumv); 
eta0=eta0';

eta1s(i,1)=eta0(1);
eta2s(i,1)=eta0(2);
eta3s(i,1)=eta0(3);
eta4s(i,1)=eta0(4);
eta5s(i,1)=eta0(5);
eta6s(i,1)=eta0(6);
eta7s(i,1)=eta0(7);
eta8s(i,1)=eta0(8);
eta9s(i,1)=eta0(9);
eta10s(i,1)=eta0(10);

%Sample sigmam
 sums=0;

for t=1:T
   
    
    fivt=fiv{t}; omegat=omega0{t};   w5gt=W5g{t};  xit=xi0{t};
     wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
   
    for ii=1:n
        for jj=1:n
        if (ii~=jj)
         r=[fivt(ii,jj),B1(ii,t),B1(jj,t),wea1t(ii,jj),wea2t(ii,jj),wea3t(ii,jj),wea4t(ii,jj),X47(ii),X47(jj),M(ii,jj)];  
       
        ss=w5gt(ii,jj)-r*eta0-omegat(ii)*kappa40-xit(jj)*zeta40;
  sums=sums+ss'*ss;   
        end
        
        

        end
    end
end

 
 
 ap=a+n*(n-1)*T/2;
 bp=b+0.5*sums;
sigmam0= 1./gamrnd(ap,1/bp);  
sigmams(i,1)=sigmam0;   


%sample kappa40

 sumk=0; sumkv=0;

for t=1:T
  
    fivt=fiv{t}; omegat=omega0{t};   w5gt=W5g{t};  xit=xi0{t};
     wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
     
    for ii=1:n
        for jj=1:n
        if (ii~=jj)
          r=[fivt(ii,jj),B1(ii,t),B1(jj,t),wea1t(ii,jj),wea2t(ii,jj),wea3t(ii,jj),wea4t(ii,jj),X47(ii),X47(jj),M(ii,jj)]; 
        yy=w5gt(ii,jj)-r*eta0-xit(jj)*zeta40;
  sumk=sumk+omegat(ii)*yy/sigmam0;   
  sumkv=sumkv+omegat(ii)^2/sigmam0;
        end
        
        

        end
    end
end

 sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;
kappa40=mvnrnd(Tr,sumv); 

%  kappa40=normt_rnd(Tr,sumv,0,9999999);

kappa4s(i,1)=kappa40;



%sample kappa50

 sumk=0; sumkv=0;

for t=1:T
    
    
    fivt=fiv{t}; omegat=omega0{t};  w5gt=W5g{t};  xit=xi0{t};
     wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
     
    for ii=1:n
        for jj=1:n
        if (ii~=jj)
       r=[fivt(jj,ii),B1(jj,t),B1(ii,t),wea1t(jj,ii),wea2t(jj,ii),wea3t(jj,ii),wea4t(jj,ii),X47(jj),X47(ii),M(jj,ii)]; 
       
        yy=w5gt(jj,ii)-r*eta0-omegat(jj)*kappa40;
  sumk=sumk+xit(ii)*yy/sigmam0;   
  sumkv=sumkv+xit(ii)^2/sigmam0;
        end
        
        

        end
    end
end

 sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;
zeta40=mvnrnd(Tr,sumv); 

%  zeta40=normt_rnd(Tr,sumv,0,9999999);

zeta4s(i,1)=zeta40;


%***************
% Counterfactual 
%***************

%store residual estimates and CF cross state flows, and vaccination.
Res=zeros(n,T); taues=zeros(n,T); hFw=zeros(n,T); Resv=zeros(n,T);
mues=cell(T,1); hW=cell(T,1); hWL=cell(T,1);
for t=1:T
 mues{t}=zeros(n,n);   
 hW{t}=zeros(n,n);
 hWL{t}=zeros(n,n);
 
end


%obtain residual estimates
if (i>nomit)
    
  for t=1:T
        
        D0{t}=diag(B1(:,t)<=gamma0(1));
        M0{t}=diag(B1(:,t)<=gamma0(2));
        YL1=D0{t}*inpop*YL(:,t); YL2=(eye(n)-D0{t})*inpop*YL(:,t);
        YL3=M0{t}*inpopr*YL(:,t); YL4=(eye(n)-M0{t})*inpopr*YL(:,t);
        
        W10=W5{t};
        ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
        
        SS=eye(n)-lambda0(1)*W10;
        Yrr=SS*Y(:,t)-rho0(1)*YL1-rho0(2)*YL2-rho0(3)*YL3-rho0(4)*YL4-mu0*W5L{t}*YL(:,t);
        Res(:,t)=Yrr-ZZ*beta0-Time10(t,1)*ones(n,1)-indiv10-omega0{t}*kappa10-xi0{t}*zeta10;
        
         ZZp=[B1(:,t),Fw19L(:,t),X1(:,:,t)];
       
        taues(:,t)=Fw(:,t)-ZZp*phi0-omega0{t}*kappa30-xi0{t}*zeta30-ones(n,1)*Time30(t)-indiv30;
        
        
             ZZv=[BBL(:,t),E1(:,t),ew(:,:,t)];
        Resv(:,t)=B1(:,t)-ZZv*delta0-indiv20-omega0{t}*kappa20-xi0{t}*zeta20;
        
        
  end
  
  
  
    %obtain the residuals of the gravity equation
    
    
    for t=1:T
        fivt=fiv{t}; omegat=omega0{t};  w5gt=W5g{t}; xit=xi0{t};
        wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
        
        for ii=1:n
            for jj=1:n
                
                
                if (ii~=jj)
                    r=[fivt(ii,jj),B1(ii,t),B1(jj,t),wea1t(ii,jj),wea2t(ii,jj),wea3t(ii,jj),wea4t(ii,jj),X47(ii),X47(jj),M(ii,jj)];
                    
                    rr=w5gt(ii,jj)-r*eta0-omegat(ii)*kappa40-xit(jj)*zeta40;
                    mues{t}(ii,jj)=rr;
                end
                
                
                
            end
        end
        
    end
  
  
  
  
  YLc=zeros(n,T);   YLc(:,1)=YL(:,1); 
  
  
if sce==2
      
         
     

for t=1:T
    
    omegat=omega0{t}; xit=xi0{t};
    D0{t}=diag(B1(:,t)<=gamma0(1));  M0{t}=diag(B1(:,t)<=gamma0(2));
    
    W10=W5{t};
    ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
    SS=eye(n)-lambda0(1)*W10;
    SSi=SS\eye(n);
    xs=(rho0(1)*D0{t}*inpop+rho0(1)*(eye(n)-D0{t})*inpop+rho0(3)*M0{t}*inpopr+rho0(3)*(eye(n)-M0{t})*inpopr+mu0*W5L{t})*YLc(:,t)+ZZ*beta0+indiv10+Time10(t,1)*ones(n,1)+omegat*kappa10+xit*zeta10+Res(:,t);
    ylc=SSi*xs;

    Yc(t,i)=sum(ylc);
    
    for ii=1:n
      statecase(ii,t,i-nomit)=ylc(ii);  
         Ymean(ii,t)=Ymean(ii,t)+ylc(ii);
       stateccase(ii,i-nomit)=stateccase(ii,i-nomit)+ylc(ii);  
    end
   
    if (t<T)
        YLc(:,t+1)=ylc;
    end
    
    
    
    
end



 elseif sce==3
     
   phin0=phi0; phin0(1)=0;
   
   

    for t=1:T
        
   omegat=omega0{t}; xit=xi0{t};
           
   hFw(:,t)=[B1(:,t),Fw19L(:,t),X1(:,:,t)]*phin0+indiv30+ones(n,1)*Time30(t)+omegat*kappa30+xit*zeta30+taues(:,t);
        
  D0{t}=diag(B1(:,t)<=gamma0(1));  M0{t}=diag(B1(:,t)<=gamma0(2));
    
    W10=W5{t};
    ZZ=[hFw(:,t), X1(:,:,t),nr(:,:,t)];
    SS=eye(n)-lambda0(1)*W10;
    SSi=SS\eye(n);
    xs=(rho0(1)*D0{t}*inpop+rho0(2)*(eye(n)-D0{t})*inpop+rho0(3)*M0{t}*inpopr+rho0(4)*(eye(n)-M0{t})*inpopr+mu0*W5L{t})*YLc(:,t)+ZZ*beta0+indiv10+Time10(t,1)*ones(n,1)+omegat*kappa10+xit*zeta10+Res(:,t);
    ylc=SSi*xs;

    Yc(t,i)=sum(ylc);
    
          
    
    for ii=1:n
      statecase(ii,t,i-nomit)=ylc(ii);  
         Ymean(ii,t)=Ymean(ii,t)+ylc(ii);
        stateccase(ii,i-nomit)=stateccase(ii,i-nomit)+ylc(ii); 
    end
    
    if (t<T)
        YLc(:,t+1)=ylc;
    end
     
    end
     
 elseif sce==4
    
     
     
     
     %predicting network flow
      etan0=eta0; etan0(2)=0; etan0(3)=0;
      
   
      
      for t=1:T
          fivt=fiv{t}; omegat=omega0{t};   xit=xi0{t};
           wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
          
           
          for ii=1:n
               for jj=1:n
               
                   
                   if (ii~=jj)
                       r=[fivt(ii,jj),B1(ii,t),B1(jj,t),wea1t(ii,jj),wea2t(ii,jj),wea3t(ii,jj),wea4t(ii,jj),X47(ii),X47(jj),M(ii,jj)];
                       
                       ww=r*etan0+omegat(ii)*kappa40+xit(jj)*zeta40+mues{t}(ii,jj);
                       hW{t}(ii,jj)=ww;
                      
                   end
                   
                   
                   
               end
          end 
           
         
      end
      
      
      for t=1:T
         if t==1
             hWL{t}=W5L{t};
         else
             hWL{t}=hW{t-1};
         end
          
      end
      
      rsum=zeros(T,1);  rsuml=zeros(T,1);

for t=1:T
    rsum(t,1)=max(sum(abs(hW{t}),2));
    rsuml(t,1)=max(sum(abs(hWL{t}),2));
end

factor=max(max(rsum),max(rsuml));
      
  for t=1:T
    %row-normalized population flow from Wuhan
 ww=hW{t};
ww=ww./factor;
   hW{t}=ww;
 end

for t=1:T

   w1=hWL{t};

w1=w1./factor;
      
    hWL{t}=w1;
end
%     
      
       
    
         
for t=1:T
           
  omegat=omega0{t}; xit=xi0{t};
           
  
        
  D0{t}=diag(B1(:,t)<=gamma0(1));  M0{t}=diag(B1(:,t)<=gamma0(2));
    
    W10=hW{t};
    ZZ=[Fw(:,t), X1(:,:,t),nr(:,:,t)];
    SS=eye(n)-lambda0(1)*W10;
    SSi=SS\eye(n);
    xs=(rho0(1)*D0{t}*inpop+rho0(2)*(eye(n)-D0{t})*inpop+rho0(3)*M0{t}*inpopr+rho0(4)*(eye(n)-M0{t})*inpopr+mu0*hWL{t})*YLc(:,t)+ZZ*beta0+Time10(t,1)*ones(n,1)+indiv10+omegat*kappa10+xit*zeta10+Res(:,t);
    ylc=SSi*xs;

    Yc(t,i)=sum(ylc);
    
          
    
  for ii=1:n
      statecase(ii,t,i-nomit)=ylc(ii);  
      stateccase(ii,i-nomit)=stateccase(ii,i-nomit)+ylc(ii);
         Ymean(ii,t)=Ymean(ii,t)+ylc(ii);
         
         
         % outflow=sum(abs(W10),1); inflow=sum(abs(W10),2);
          
           outflow=sum(W10,1); inflow=sum(W10,2);
         
     stateoutflow(ii,t,i-nomit)=outflow(ii);    
     outflowm(ii,t)=outflowm(ii,t)+outflow(ii);
     
        stateinflow(ii,t,i-nomit)=inflow(ii);      
          inflowm(ii,t)=inflowm(ii,t)+inflow(ii);
        
  end
    
    if (t<T)
        YLc(:,t+1)=ylc;
    end
          
end
    
 elseif sce==5
     
     
     
     %predicting network flow
      etan0=eta0; etan0(2)=0; etan0(3)=0;
      
   
      
      for t=1:T
          fivt=fiv{t}; omegat=omega0{t};   xit=xi0{t};
           wea1t=wea1{t}; wea2t=wea2{t}; wea3t=wea3{t}; wea4t=wea4{t};
           
          for ii=1:n
               for jj=1:n
               
                   
                   if (ii~=jj)
                       r=[fivt(ii,jj),B1(ii,t),B1(jj,t),wea1t(ii,jj),wea2t(ii,jj),wea3t(ii,jj),wea4t(ii,jj),X47(ii),X47(jj),M(ii,jj)];
                       
                       ww=r*etan0+omegat(ii)*kappa40+xit(jj)*zeta40+mues{t}(ii,jj);
                       hW{t}(ii,jj)=ww;
                      
                   end
                   
                   
                   
               end
           end 
           
          
      end
      
      
      for t=1:T
         if t==1
             hWL{t}=W5L{t};
         else
             hWL{t}=hW{t-1};
         end
          
      end
      
      rsum=zeros(T,1);  rsuml=zeros(T,1);

for t=1:T
    rsum(t,1)=max(sum(abs(hW{t}),2));
    rsuml(t,1)=max(sum(abs(hWL{t}),2));
end

factor=max(max(rsum),max(rsuml));
      
  for t=1:T
    %row-normalized population flow from Wuhan
 ww=hW{t};
ww=ww./factor;
   hW{t}=ww;
   Wp=Wp+hW{t};
 end

for t=1:T

   w1=hWL{t};

w1=w1./factor;
      
    hWL{t}=w1;
end 
  

phin0=phi0; phin0(1)=0;
       
      
         
       for t=1:T
%          ZZp=[B1(:,t),Fw19L(:,t),X1(:,:,t)];
%        
%         taues(:,t)=Fw(:,t)-ZZp*phi0-omega0{t}*kappa30-xi0{t}*zeta30-ones(n,1)*Time30(t)-indiv30;
%            
    omegat=omega0{t}; xit=xi0{t};
           
    hFw(:,t)=[B1(:,t),Fw19L(:,t),X1(:,:,t)]*phin0+indiv30+ones(n,1)*Time30(t)+omegat*kappa30+xit*zeta30+taues(:,t);
        
  D0{t}=diag(B1(:,t)<=gamma0(1));  M0{t}=diag(B1(:,t)<=gamma0(2));
    
    W10=hW{t};
    ZZ=[hFw(:,t), X1(:,:,t),nr(:,:,t)];
    SS=eye(n)-lambda0(1)*W10;
    SSi=SS\eye(n);
    xs=(rho0(1)*D0{t}*inpop+rho0(1)*(eye(n)-D0{t})*inpop+rho0(3)*M0{t}*inpopr+rho0(3)*(eye(n)-M0{t})*inpopr+mu0*hWL{t})*YLc(:,t)+ZZ*beta0+Time10(t,1)*ones(n,1)+indiv10+omegat*kappa10+xit*zeta10+Res(:,t);
    ylc=SSi*xs;

    Yc(t,i)=sum(ylc);
    
          
    
    for ii=1:n
      statecase(ii,t,i-nomit)=ylc(ii);  
      stateccase(ii,i-nomit)=stateccase(ii,i-nomit)+ylc(ii);
         Ymean(ii,t)=Ymean(ii,t)+ylc(ii);
         
         
%           outflow=sum(abs(W10),1); inflow=sum(abs(W10),2);
          
          outflow=sum(W10,1); inflow=sum(W10,2);  
         
     stateoutflow(ii,t,i-nomit)=outflow(ii);    
     outflowm(ii,t)=outflowm(ii,t)+outflow(ii);
     
        stateinflow(ii,t,i-nomit)=inflow(ii);      
          inflowm(ii,t)=inflowm(ii,t)+inflow(ii);
        
    end
    
     if (t<T)
        YLc(:,t+1)=ylc;
    end    
          
       end
    
 
end 

end  
      
      
      
      
      


end



%compute the UB and LB of CF new case, CF vax for each state
for t=1:T
   for ii=1:n
      aa=zeros(nit-nomit,1);
      vv=zeros(nit-nomit,1);
        in=zeros(nit-nomit,1);
      out=zeros(nit-nomit,1);
      for iii=1:nit-nomit
        aa(iii,1)=statecase(ii,t,iii);
        vv(iii,1)=statevax(ii,t,iii);
          in(iii,1)=stateoutflow(ii,t,iii);
        out(iii,1)=stateinflow(ii,t,iii);
      end
      
%       aa=aa';
      bound=hpdi(aa,0.95);
      Yub(ii,t)=bound(1);
       Ylb(ii,t)=bound(2);
       
        boundv=hpdi(vv,0.95);
      Vub(ii,t)=boundv(1);
       Vlb(ii,t)=boundv(2);
       
        boundi=hpdi(in,0.95);
       outflowu(ii,t)=boundi(1);
       outflowl(ii,t)=boundi(2);
       
        boundo=hpdi(out,0.95);
       inflowu(ii,t)=boundo(1);
       inflowl(ii,t)=boundo(2);
       
       
   end
end  


% compute the CF cumulative at the end of the sample for each state
for ii=1:n
 cc=stateccase(ii,:);
    cc=cc';
    Ycmean(ii,1)=mean(cc);
       bound=hpdi(cc,0.95);
      Ycub(ii,1)=bound(1);
       Yclb(ii,1)=bound(2);
    
end



    
 



results.gamma1=gamma1s(nomit+1:nit,1);
results.gamma2=gamma2s(nomit+1:nit,1);
 results.lambda=lambdas(nomit+1:nit,1);
 results.rho1=rho1s(nomit+1:nit,1);
  results.rho2=rho2s(nomit+1:nit,1);
  results.rho3=rho3s(nomit+1:nit,1);
  results.rho4=rho4s(nomit+1:nit,1);
  results.rhod=rhods(nomit+1:nit,1);
  results.mu=mus(nomit+1:nit,1);
  results.beta1=beta1s(nomit+1:nit,1);
  results.beta2=beta2s(nomit+1:nit,1);
  results.beta3=beta3s(nomit+1:nit,1);
results.beta4=beta4s(nomit+1:nit,1);
  results.beta5=beta5s(nomit+1:nit,1);
  results.beta6=beta6s(nomit+1:nit,1);
   results.betap=betaps(nomit+1:nit,1);
    results.eta1=eta1s(nomit+1:nit,1);
      results.eta2=eta2s(nomit+1:nit,1);
       results.eta3=eta3s(nomit+1:nit,1);
         results.eta4=eta4s(nomit+1:nit,1);
            results.eta5=eta5s(nomit+1:nit,1);
        results.eta6=eta6s(nomit+1:nit,1);
      results.eta7=eta7s(nomit+1:nit,1);
 results.eta8=eta8s(nomit+1:nit,1);
 results.eta9=eta9s(nomit+1:nit,1);


  results.delta1=delta1s(nomit+1:nit,1);
  results.delta2=delta2s(nomit+1:nit,1);
  
  results.delta3=delta3s(nomit+1:nit,1);
  results.delta4=delta4s(nomit+1:nit,1);
  
  results.delta5=delta5s(nomit+1:nit,1);
    results.delta6=delta6s(nomit+1:nit,1);
 
   results.phi1=phi1s(nomit+1:nit,1);
    results.phi2=phi2s(nomit+1:nit,1);
       results.phi3=phi3s(nomit+1:nit,1);
    results.phi4=phi4s(nomit+1:nit,1);
    results.phi5=phi5s(nomit+1:nit,1);
      results.phi6=phi6s(nomit+1:nit,1);
    results.phip=phips(nomit+1:nit,1);

 results.kappa1=kappa1s(nomit+1:nit,1);
 results.kappa2=kappa2s(nomit+1:nit,1);
 results.kappa3=kappa3s(nomit+1:nit,1);
 results.sigmav=sigmavs(nomit+1:nit,1);
 results.sigmae=sigmaes(nomit+1:nit,1);
  results.sigmat=sigmats(nomit+1:nit,1);
results.sigmam=sigmams(nomit+1:nit,1);

results.Yc=Yc(:,nomit+1:nit);
results.Vc=Vc(:,nomit+1:nit);
results.Ymean=Ymean./(nit-nomit);
results.Yub=Yub;
results.Ylb=Ylb;
results.Ycmean=Ycmean;
results.Ycub=Ycub;
results.Yclb=Yclb;
results.Vmean=Vmean./(nit-nomit);
results.Vub=Vub;
results.Vlb=Vlb;
results.Wp=Wp/(nit-nomit);
results.outflowm=outflowm/(nit-nomit);
results.outflowu=outflowu;
results.outflowl=outflowl;
results.inflowm=inflowm/(nit-nomit);
results.inflowu=inflowu;
results.inflowl=inflowl;
%results.Yp=Yp./(nit-nomit);
%results.Ypcum=Ypcum(:,nomit+1:nit);

end

% =========================================================================
% support functions below
% =========================================================================


function [fvalue ] = likn(Yq, Ylq, Y10, Y20,Y30,Y40, Y11, Y21,Y31,Y41,nrq, X1q, fq, W10,WW,omegaq,xiq, indiv10, Time1q,  lambda0,rho0,mu0,lambda1,rho1,mu1,beta0,kappa10,zeta10,sigmav0)

betar0=beta0(2:5); 



n=length(Yq);

S0=eye(n)-lambda0*W10;
S1=eye(n)-lambda1*W10;

C0=S0*Yq-rho0(1)*Y10-rho0(2)*Y20-rho0(3)*Y30-rho0(4)*Y40-mu0*WW*Ylq-fq*beta0(1)-X1q*betar0-nrq*beta0(6)-indiv10-ones(n,1)*Time1q-omegaq*kappa10-xiq*zeta10;
C1=S1*Yq-rho1(1)*Y11-rho1(2)*Y21-rho1(3)*Y31-rho1(4)*Y41-mu1*WW*Ylq-fq*beta0(1)-X1q*betar0-nrq*beta0(6)-indiv10-ones(n,1)*Time1q-omegaq*kappa10-xiq*zeta10;

CC0=C0'*C0/(2*sigmav0);
CC1=C1'*C1/(2*sigmav0);

fvalue=(det(S1)/det(S0))*exp(-CC1+CC0);


end
