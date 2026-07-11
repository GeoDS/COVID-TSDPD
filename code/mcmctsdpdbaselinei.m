function [results] =mcmctsdpdbaselinei(nit,br, n,T,Y,YL,X1,B1,Bd1,BdL,BBdL,Xv,X31,Xstate,Xtr,M,W5,inpop,W5g,W19g,wea1,wea2,Time10,Time30,lambda0, rho0,gamma0,lgamma,ugamma,beta0, delta0,eta0,kappa0,zeta0,sigmav0,sigmae0,sigmac0,rs)
% =========================================================================
% This function conducts the MCMC estimation of a higher-order SAR model with spatial errors.
%==========================================================================

nomit=nit*br; %burn-in sample


% 
%  Brn=zeros(n,4,T);  
%store all parameters
%parameter in the SDPD equation
lambdas=zeros(nit,1);  gamma1s=zeros(nit,1); gamma2s=zeros(nit,1); gamma3s=zeros(nit,1); gamma4s=zeros(nit,1);
rho1s=zeros(nit,1); rho2s=zeros(nit,1); rho3s=zeros(nit,1); rho4s=zeros(nit,1); 
rho5s=zeros(nit,1); rho6s=zeros(nit,1); 

beta1s=zeros(nit,1); beta2s=zeros(nit,1); beta3s=zeros(nit,1); beta4s=zeros(nit,1); beta5s=zeros(nit,1);
beta6s=zeros(nit,1); sigmavs=zeros(nit,1); kappa1s=zeros(nit,1); zeta1s=zeros(nit,1);   
%  Time1s=zeros(T,nit); 

%parameter for the linear equation on vaccination
delta1s=zeros(nit,1);  delta2s=zeros(nit,1); delta3s=zeros(nit,1); delta4s=zeros(nit,1); delta5s=zeros(nit,1);
sigmaes=zeros(nit,1);  kappa2s=zeros(nit,1);   zeta2s=zeros(nit,1);   delta6s=zeros(nit,1);
% Time2s=zeros(T,nit); 


  
%parameter for gravity equation
kappa3s=zeros(nit,1);  zeta3s=zeros(nit,1); 
eta1s=zeros(nit,1); eta2s=zeros(nit,1); eta3s=zeros(nit,1); eta4s=zeros(nit,1); eta5s=zeros(nit,1);
eta6s=zeros(nit,1); eta7s=zeros(nit,1); eta8s=zeros(nit,1); eta9s=zeros(nit,1);eta10s=zeros(nit,1); 
sigmacs=zeros(nit,1); 
% Time4s=zeros(T,nit); 

%store transform data
Yr=zeros(n,T); DL0=cell(T,1); ML0=cell(T,1);  DL1=cell(T,1); ML1=cell(T,1);
DH0=cell(T,1); MH0=cell(T,1);  DH1=cell(T,1); MH1=cell(T,1);


%   
%=======================================
%prepare for the AM algorithm for lambda 
%=======================================
KKt=7;   st1=0.1/sqrt(KKt); ratio=0.05;
sum1=zeros(7,1);  sum2=zeros(7,7);


KKg=4; sg1=0.1/sqrt(KKg); sumg1=zeros(4,1); sumg2=zeros(4,4);


% % prepare to calculate the mcmc draws of direct effects
% directs=zeros(nit,1);
% 
% Ynf=sum(sum(Y)); 
% Yseed=zeros(T,1);
% 
% for t=1:T
%    Yseed(t,1)=sum(Y(index,t));
% end
% 
% Yseedf=sum(Yseed);




%==========================
% value of prior parameters
%===========================
TT=200; %initial period of AM algorithm 
P=1000000;  % for beta1 and beta2   
a=0.001;b=0.001;  %for sigma
ae=0.001; be=0.001;
ac=0.001; bc=0.001;
theta0=[lambda0;rho0];


inpopr=eye(n)-inpop;

omega0=cell(T,1); xi0=cell(T,1);
for t=1:T
   omega0{t}=zeros(n,1); 
    xi0{t}=zeros(n,1);
end

kappa10=kappa0(1); kappa20=kappa0(2); kappa30=kappa0(3); 
zeta10=zeta0(1); zeta20=zeta0(2); zeta30=zeta0(3); 


for t=1:T
      DL0{t}=diag(B1(:,t)<=gamma0(1));
      ML0{t}=diag(B1(:,t)>gamma0(2));
      DH0{t}=diag(B1(:,t)<=gamma0(3));
      MH0{t}=diag(B1(:,t)>gamma0(4));
      
  YL1=DL0{t}*inpop*YL(:,t); YL2=(eye(n)-DL0{t}-ML0{t})*inpop*YL(:,t);
  YL3=ML0{t}*inpop*YL(:,t);
  
    YL4=DH0{t}*inpopr*YL(:,t); YL5=(eye(n)-DH0{t}-MH0{t})*inpopr*YL(:,t);
  YL6=MH0{t}*inpopr*YL(:,t);

    
  W10=W5{t}; 
   

    SS=eye(n)-lambda0*W10;
 Yr(:,t)=SS*Y(:,t)-rho0(1)*YL1-rho0(2)*YL2-rho0(3)*YL3-rho0(4)*YL4-rho0(5)*YL5-rho0(6)*YL6; 
end


    
for i=1: nit
    
     if mod(i,20)==0
           fprintf('%d\n',i)
          fprintf('lambda0:   %.3f %.3f %.3f\n', lambda0);
          fprintf('rho10:   %.3f %.3f %.3f\n', rho0(1));
          fprintf('rho20:   %.3f %.3f %.3f\n', rho0(2));
          fprintf('rho30:   %.3f %.3f %.3f\n', rho0(3));
          fprintf('rho40:   %.3f %.3f %.3f\n', rho0(4));
             fprintf('rho50:   %.3f %.3f %.3f\n', rho0(5));
          fprintf('rho60:   %.3f %.3f %.3f\n', rho0(6));
          fprintf('gamma10: %.3f %.3f %.3f\n', gamma0(1));
             fprintf('gamma20: %.3f %.3f %.3f\n', gamma0(2));
              fprintf('gamma30: %.3f %.3f %.3f\n', gamma0(3));
             fprintf('gamma40: %.3f %.3f %.3f\n', gamma0(4));
     

%              
%                fprintf('beta10: %.3f %.3f %.3f\n', beta0(1));
          

           fprintf('delta10: %.3f %.3f %.3f\n', delta0(1));
              fprintf('delta20: %.3f %.3f %.3f\n', delta0(2));
            
%                 
%             fprintf('eta10: %.3f %.3f %.3f\n', eta0(1));
%              fprintf('eta20: %.3f %.3f %.3f\n', eta0(2));
%               fprintf('eta30: %.3f %.3f %.3f\n', eta0(3));
%               fprintf('sigmac0: %.3f %.3f %.3f\n', sigmac0);


     end
     
     
     % sample destination latent omega
N=n-1;
sigma0i=[1/sigmav0, 0, zeros(1,N) ; 0, 1/sigmae0, zeros(1,N) ;  zeros(N,1),zeros(N,1),eye(N)/sigmac0]; kappar0=[kappa10; kappa20;ones(N,1)*kappa30];
 %kappab=kron(kappar0,eye(n)); 
  kappab=kron(eye(n),kappar0);
 %Sigma0i=kron(sigma0i,eye(n));
 Sigma0i=kron(eye(n),sigma0i);
 Sigmaw=(eye(n)+kappab'*Sigma0i*kappab)\eye(n);
 
for t=1:T
  YY=zeros((N+2)*n,1); fivt=W19g{t}; w5gt=W5g{t}; xit=xi0{t};
  wea1t=wea1{t}; wea2t=wea2{t};

 
  
   for ii=1:n
    yy3=[]; 
     zz1=[X1(ii,:,t),Xstate(ii,:),Xtr(ii,:,t)];
   
     zz2=[BBdL(ii,t),BdL(ii,t),Xv(ii,:,t),Xstate(ii,:)];
     
     
     
     yy1=Yr(ii,t)-zz1*beta0-Time10(t)-xit(ii)*zeta10;
     yy2=Bd1(ii,t)-zz2*delta0-xit(ii)*zeta20;


     for jj=1:n
        if (ii~=jj)
          r=[fivt(ii,jj),B1(ii,t),B1(jj,t),X31(ii),X31(jj),M(ii,jj),wea1t(ii,jj),wea2t(ii,jj)]; 
        gg=w5gt(ii,jj)-r*eta0-xit(jj)*zeta30-Time30(t);
         yy3=[yy3;gg];
        end
     end
     % Xstate=[X31,X42,X45,log(X46),X47/100,X48/100,X49/100,X53/100,X54/100,X55/100,X56/100];
     




     
     yy=[yy1;yy2;yy3];
     YY((ii-1)*50+1:ii*50)=yy;

   end

   Tw=Sigmaw*kappab'*(Sigma0i)*YY;
   omegaft=mvnrnd(Tw,Sigmaw);
   omegafn=omegaft';

   omega0{t}=omegafn;
 
end







     % sample origin latent xi
N=n-1;
sigma0i=[1/sigmav0, 0, zeros(1,N) ; 0, 1/sigmae0,zeros(1,N) ;  zeros(N,1),zeros(N,1),eye(N)/sigmac0]; zetar0=[zeta10; zeta20;ones(N,1)*zeta30];
 %kappab=kron(kappar0,eye(n)); 
  zetab=kron(eye(n),zetar0);
 %Sigma0i=kron(sigma0i,eye(n));
 Sigma0i=kron(eye(n),sigma0i);
 Sigmax=(eye(n)+zetab'*Sigma0i*zetab)\eye(n);
 
for t=1:T
  YY=zeros((N+2)*n,1); fivt=W19g{t}; w5gt=W5g{t}; omegat=omega0{t};
  wea1t=wea1{t}; wea2t=wea2{t}; 
  
   for ii=1:n
    yy3=[]; 
     zz1=[X1(ii,:,t),Xstate(ii,:),Xtr(ii,:,t)];
   
     zz2=[BBdL(ii,t),BdL(ii,t),Xv(ii,:,t),Xstate(ii,:)];
     
    
    
     yy1=Yr(ii,t)-zz1*beta0-Time10(t)-omegat(ii)*kappa10;
     yy2=Bd1(ii,t)-zz2*delta0-omegat(ii)*kappa20;
    

     for jj=1:n
        if (ii~=jj)
       r=[fivt(jj,ii),B1(jj,t),B1(ii,t),X31(jj),X31(ii),M(jj,ii),wea1t(jj,ii),wea2t(jj,ii)]; 
        gg=w5gt(jj,ii)-r*eta0-omegat(jj)*kappa30-Time30(t);
         yy3=[yy3;gg];
        end
     end
     
   




     
     yy=[yy1;yy2;yy3];
     YY((ii-1)*50+1:ii*50)=yy;

   end

   Tx=Sigmax*zetab'*(Sigma0i)*YY;
   xift=mvnrnd(Tx,Sigmax);
   xifn=xift';

   xi0{t}=xifn;
 
end



%======================================
% sample the SDPD model for infection
%=====================================
     
    
  if (i<=2*TT)
      accept=0;
      theta1=mvnrnd(theta0,st1^2*eye(7));
      theta1=theta1';
   lambda1=theta1(1);
      rho1=theta1(2:7);
      
  
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));
     
     
     
      
      while (accept==0) %reject bounds on lambda1 and rho1
          if (lambdam*rs<1 ) && (lambdam*rs+rhom<1)
              accept=1;
          else
           theta1=mvnrnd(theta0,st1^2*eye(7));
      theta1=theta1';
   lambda1=theta1(1);
      rho1=theta1(2:7);
  

     
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));

     
          end
      end
  end
  
  if (i>2*TT)
      accept=0;
      if (i<=nomit)
      vvarr=(1-ratio)^2*2.38^2*varr/KKt+ratio^2*st1^2*eye(7);
      end
      theta1=mvnrnd(theta0,vvarr);
      theta1=theta1';
 
  
   lambda1=theta1(1);
      rho1=theta1(2:7);

     
     
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));

      
      while (accept==0) %reject bounds on lambda1
          if (lambdam*rs<1 ) && (lambdam*rs+rhom<1)
              accept=1;
          else
      theta1=mvnrnd(theta0,vvarr);
    theta1=theta1';
 
    lambda1=theta1(1);
      rho1=theta1(2:7);

     
     
      
      lambdam=max(abs(lambda1));
      rhom=max(abs(rho1));

     
          end
      end
  end

   
 
 frpp=1;
 

  
  for t=1: T
    
       DL0{t}=diag(B1(:,t)<=gamma0(1));
      ML0{t}=diag(B1(:,t)>gamma0(2));
      
       DH0{t}=diag(B1(:,t)<=gamma0(3));
      MH0{t}=diag(B1(:,t)>gamma0(4));
      
  YL10=DL0{t}*inpop*YL(:,t); YL20=(eye(n)-DL0{t}-ML0{t})*inpop*YL(:,t);
  YL30=ML0{t}*inpop*YL(:,t);
  YL40=DH0{t}*inpopr*YL(:,t); YL50=(eye(n)-DH0{t}-MH0{t})*inpopr*YL(:,t);
  YL60=MH0{t}*inpopr*YL(:,t);
  
  W10=W5{t}; 
   
   
 frpp=likn(Y(:,t), YL10,YL20,YL30,YL40,YL50,YL60, YL10,YL20, YL30, YL40,YL50,YL60, X1(:,:,t), Xstate,Xtr(:,:,t),W10,omega0{t},xi0{t},Time10(t),lambda0,rho0,lambda1,rho1,beta0,kappa10,zeta10,sigmav0)*frpp;
  end
    

    
  % Determine the transition probability
    Acceptr=min(1,frpp);
   
    % Draw from uniform(0,1)
    u1=rand(1,1); 
    
   % Transition to candidate delta11 with probability acceptr
     if (Acceptr>u1) 
        lambda_1=lambda1;
        rho_1=rho1;
 
 
       
     else
       lambda_1=lambda0;
        rho_1=rho0;

  
       
         
     end
    
    % Store the value and continue
    rho0=rho_1;
  rho1s(i,1)=rho0(1);
  rho2s(i,1)=rho0(2);
  rho3s(i,1)=rho0(3);
  rho4s(i,1)=rho0(4);
   rho5s(i,1)=rho0(5);
  rho6s(i,1)=rho0(6);
       

     
  
     lambda0=lambda_1;
     lambdas(i,1)=lambda0;
    

      

      
  %calculate the empirical covariance
    theta0=[lambda0;rho0];
    
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
        gamma1=mvnrnd(gamma0,sg1^2*eye(4));
        gamma1=gamma1';
     
        while (accept==0) %reject bounds on lambda1 and rho1
            if  (gamma1(1)>lgamma)&& (gamma1(2)<ugamma ) && (gamma1(2)>gamma1(1))&& (gamma1(3)>lgamma)&&(gamma1(4)<ugamma)&& (gamma1(4)>gamma1(3))
                accept=1;
            else
                gamma1=mvnrnd(gamma0,sg1^2*eye(4));
                 gamma1=gamma1';
            end
        end
    end
    
    if (i>2*TT)
        accept=0;
        if (i<=nomit)
            vvarrg=(1-ratio)^2*2.38^2*varrg/KKg+ratio^2*sg1^2;
        end
        gamma1=mvnrnd(gamma0,vvarrg);
         gamma1=gamma1';
        
        while (accept==0) %reject bounds on lambda1
            if  (gamma1(1)>lgamma)&& (gamma1(2)<ugamma ) && (gamma1(2)>gamma1(1))&& (gamma1(3)>lgamma)&&(gamma1(4)<ugamma)&& (gamma1(4)>gamma1(3))
               % 
               %   (gamma1(1)<gamma1(2) ) && (gamma1(1)>lgamma) && (gamma1(2)<ugamma )
                accept=1;
            else
                gamma1=mvnrnd(gamma0,vvarrg);
                 gamma1=gamma1';
            end
        end
    end
    
    
    frpp=1;
    
    for t=1: T
       DL0{t}=diag(B1(:,t)<=gamma0(1));
        DL1{t}=diag(B1(:,t)<=gamma1(1));
        
        ML0{t}=diag(B1(:,t)>gamma0(2));
        ML1{t}=diag(B1(:,t)>gamma1(2));
        
         DH0{t}=diag(B1(:,t)<=gamma0(3));
        DH1{t}=diag(B1(:,t)<=gamma1(3));
        
        MH0{t}=diag(B1(:,t)>gamma0(4));
        MH1{t}=diag(B1(:,t)>gamma1(4));
        
        
        W10=W5{t}; 
        
        YL10=DL0{t}*inpop*YL(:,t); YL20=(eye(n)-DL0{t}-ML0{t})*inpop*YL(:,t);
        YL30=ML0{t}*inpop*YL(:,t); 
        
        
        YL40=DH0{t}*inpopr*YL(:,t);    YL50=(eye(n)-DH0{t}-MH0{t})*inpopr*YL(:,t);   
       YL60=MH0{t}*inpopr*YL(:,t); 
        
      YL11=DL1{t}*inpop*YL(:,t); YL21=(eye(n)-DL1{t}-ML1{t})*inpop*YL(:,t);
        YL31=ML1{t}*inpop*YL(:,t); 
        
        
        YL41=DH1{t}*inpopr*YL(:,t);    YL51=(eye(n)-DH1{t}-MH1{t})*inpopr*YL(:,t);   
       YL61=MH1{t}*inpopr*YL(:,t); 
        
   
        
        frpp=likn(Y(:,t), YL10,YL20,YL30,YL40,YL50,YL60,YL11,YL21,YL31,YL41,YL51,YL61,X1(:,:,t),Xstate,Xtr(:,:,t),W10,omega0{t},xi0{t},Time10(t),lambda0,rho0,lambda0,rho0,beta0,kappa10,zeta10,sigmav0)*frpp;
    end
    
%   frpp=frpp*(ugamma-gamma0(2))/(ugamma-gamma1(2));
     frpp=frpp*(ugamma-gamma0(1))/(ugamma-gamma1(1))*((ugamma-gamma0(3))/(ugamma-gamma1(3)));
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
     gamma3s(i,1)=gamma0(3);
    gamma4s(i,1)=gamma0(4);
    
    
    
    
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
    DL0{t}=diag(B1(:,t)<=gamma0(1));
      ML0{t}=diag(B1(:,t)>gamma0(2));
      DH0{t}=diag(B1(:,t)<=gamma0(3));
      MH0{t}=diag(B1(:,t)>gamma0(4));
      
  YL1=DL0{t}*inpop*YL(:,t); YL2=(eye(n)-DL0{t}-ML0{t})*inpop*YL(:,t);
  YL3=ML0{t}*inpop*YL(:,t);
  
    YL4=DH0{t}*inpopr*YL(:,t); YL5=(eye(n)-DH0{t}-MH0{t})*inpopr*YL(:,t);
  YL6=MH0{t}*inpopr*YL(:,t);
    
   

    SS=eye(n)-lambda0(1)*W5{t};
 Yr(:,t)=SS*Y(:,t)-rho0(1)*YL1-rho0(2)*YL2-rho0(3)*YL3-rho0(4)*YL4-rho0(5)*YL5-rho0(6)*YL6; 
    
end
    

   

    

   
    

  %sample beta0=[delta0; beta10; beta20]
 sumb=zeros(111,1); sumbv=zeros(111,111);
 
for t=1:T
 ZZ=[X1(:,:,t),Xstate,Xtr(:,:,t)];
 sumb=sumb+ZZ'*(Yr(:,t)-ones(n,1)*Time10(t)-omega0{t}*kappa10-xi0{t}*zeta10)/sigmav0;
 sumbv=sumbv+ZZ'*ZZ/sigmav0;

end

sumv=(sumbv+eye(111)/(P))\eye(111);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumb;

betat=mvnrnd(Tr,sumv); 


beta0=betat';
beta1s(i,1)=beta0(1);
beta2s(i,1)=beta0(2);
beta3s(i,1)=beta0(3);
beta4s(i,1)=beta0(4);
beta5s(i,1)=beta0(5);
beta6s(i,1)=beta0(6);  


 

%Sample kappa10
sumk=0; sumkv=0;
for t=1:T
 ZZ=[X1(:,:,t),Xstate,Xtr(:,:,t)];
  omegat=omega0{t};  xit=xi0{t};
sumk=sumk+omegat'*(Yr(:,t)-ones(n,1)*Time10(t)-ZZ*beta0-xit*zeta10)/sigmav0;
 sumkv=sumkv+omegat'*omegat/sigmav0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;

kappa10=mvnrnd(Tr,sumv); 

kappa1s(i,1)=kappa10;

%Sample zeta10
sumk=0; sumkv=0;
for t=1:T
 ZZ=[X1(:,:,t),Xstate,Xtr(:,:,t)];
  omegat=omega0{t};  xit=xi0{t};
sumk=sumk+xit'*(Yr(:,t)-ones(n,1)*Time10(t)-ZZ*beta0-omegat*kappa10)/sigmav0;
 sumkv=sumkv+xit'*xit/sigmav0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;

zeta10=mvnrnd(Tr,sumv); 

zeta1s(i,1)=zeta10;




%Sample sigmav
 sums=0;
 for t=1:T
  ZZ=[X1(:,:,t),Xstate,Xtr(:,:,t)];
    omegat=omega0{t}; xit=xi0{t};
  ss=Yr(:,t)-ZZ*beta0-ones(n,1)*Time10(t)-omegat*kappa10-xit*zeta10;
  sums=sums+ss'*ss;   
 end
 
 
 ap=a+n*T/2;
 bp=b+0.5*sums;
sigmav0= 1./gamrnd(ap,1/bp);  
sigmavs(i,1)=sigmav0;   


 %Gibbs sampling for time-fixed effects
    
    for t=2:T
      ZZ=[X1(:,:,t),Xstate,Xtr(:,:,t)];
          omegat=omega0{t}; xit=xi0{t};
        Yf=Yr(:,t)-ZZ*beta0-omegat*kappa10-xit*zeta10;
        
        Sigmaf=(0.001+n/(sigmav0))^(-1);
       
        Tf=Sigmaf*ones(n,1)'*Yf/(sigmav0);
        
        time0=mvnrnd(Tf,Sigmaf);
        
        Time10(t,1)=time0;
       
    end
    
    
  


%*****************************
% Sample parameters in the linear panel on vaccination
%******************************
  sume=zeros(15,1); sumev=zeros(15,15);
 
for t=1:T
 ZZ=[BBdL(:,t),BdL(:,t),Xv(:,:,t),Xstate];
    omegat=omega0{t}; xit=xi0{t};
 sume=sume+ZZ'*(Bd1(:,t)-omegat*kappa20-xit*zeta20)/sigmae0;
 sumev=sumev+ZZ'*ZZ/sigmae0;

end

sumv=(sumev+eye(15)/(P))\eye(15);
sumv=mean(cat(3,sumv,sumv'),3); 

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
   ZZ=[BBdL(:,t),BdL(:,t),Xv(:,:,t),Xstate];
     omegat=omega0{t}; xit=xi0{t}; 
 ss=Bd1(:,t)-ZZ*delta0-omegat*kappa20-xit*zeta20;
  sums=sums+ss'*ss;   
 end
 
 
 aep=ae+n*T/2;
 bep=be+0.5*sums;
sigmae0= 1./gamrnd(aep,1/bep);  
sigmaes(i,1)=sigmae0;   

 


 
% sample kappa20
sumk=0; sumkv=0;
for t=1:T
   ZZ=[BBdL(:,t),BdL(:,t),Xv(:,:,t),Xstate];
  omegat=omega0{t}; xit=xi0{t};
  sumk=sumk+omegat'*(Bd1(:,t)-ZZ*delta0-xit*zeta20)/sigmae0;
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
      ZZ=[BBdL(:,t),BdL(:,t),Xv(:,:,t),Xstate];
  omegat=omega0{t}; xit=xi0{t};
  sumk=sumk+xit'*(Bd1(:,t)-ZZ*delta0-omegat*kappa20)/sigmae0;
 sumkv=sumkv+xit'*xit/sigmae0; 
    
    
end

sumv=(sumkv+eye(1)/(P))\eye(1);
%sumv=mean(cat(3,sumv,sumv'),3); 

Tr=sumv*sumk;
zeta20=mvnrnd(Tr,sumv); 
zeta2s(i,1)=zeta20;










  

%**************************
% gravity equation
%*************************

sumeta=zeros(8,1); sumetav=zeros(8,8);
for t=1:T

    fivt=W19g{t}; omegat=omega0{t};  w5gt=W5g{t}; xit=xi0{t};
    wea1t=wea1{t}; wea2t=wea2{t};  
  
    for ii=1:n
        for jj=1:n
     
            
        if (ii~=jj)
         r=[fivt(ii,jj),B1(ii,t),B1(jj,t),X31(ii),X31(jj),M(ii,jj),wea1t(ii,jj),wea2t(ii,jj)]; 
        
         
         
         sumeta=sumeta+r'*(w5gt(ii,jj)-omegat(ii)*kappa30-xit(jj)*zeta30-Time30(t))/sigmac0;
        sumetav=sumetav+r'*r/sigmac0;
        
        end
        
        

        end
    end
end

sumv=(sumetav+eye(8)/(P))\eye(8);

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



%Sample sigmac
 sums=0;

for t=1:T
   
    
    fivt=W19g{t}; omegat=omega0{t};   w5gt=W5g{t};  xit=xi0{t};
     wea1t=wea1{t}; wea2t=wea2{t}; 
   
    for ii=1:n
        for jj=1:n
        if (ii~=jj)
          r=[fivt(ii,jj),B1(ii,t),B1(jj,t),X31(ii),X31(jj),M(ii,jj),wea1t(ii,jj),wea2t(ii,jj)]; 
       
        ss=w5gt(ii,jj)-r*eta0-omegat(ii)*kappa30-xit(jj)*zeta30-Time30(t);
  sums=sums+ss'*ss;   
        end
        
        

        end
    end
end

 
 acp=ac+n*(n-1)*T/2;
 bcp=bc+0.5*sums;
sigmac0= 1./gamrnd(acp,1/bcp);  
sigmacs(i,1)=sigmac0;   






%sample kappa30

 sumk=0; sumkv=0;

for t=1:T
  
    fivt=W19g{t}; omegat=omega0{t};   w5gt=W5g{t};  xit=xi0{t};
     wea1t=wea1{t}; wea2t=wea2{t};  
     
    for ii=1:n
        for jj=1:n
        if (ii~=jj)
          r=[fivt(ii,jj),B1(ii,t),B1(jj,t),X31(ii),X31(jj),M(ii,jj),wea1t(ii,jj),wea2t(ii,jj)]; 
        yy=w5gt(ii,jj)-r*eta0-xit(jj)*zeta30-Time30(t);
  sumk=sumk+omegat(ii)*yy/sigmac0;   
  sumkv=sumkv+omegat(ii)^2/sigmac0;
        end
        
        

        end
    end
end

 sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;
kappa30=mvnrnd(Tr,sumv); 

%  kappa30=normt_rnd(Tr,sumv,0,9999999);

kappa3s(i,1)=kappa30;



%sample zeta30

 sumk=0; sumkv=0;

for t=1:T
    
    
    fivt=W19g{t}; omegat=omega0{t};  w5gt=W5g{t};  xit=xi0{t};
     wea1t=wea1{t}; wea2t=wea2{t}; 
     
    for ii=1:n
        for jj=1:n
        if (ii~=jj)
      r=[fivt(jj,ii),B1(jj,t),B1(ii,t),X31(jj),X31(ii),M(jj,ii),wea1t(jj,ii),wea2t(jj,ii)]; 
       
        yy=w5gt(jj,ii)-r*eta0-omegat(jj)*kappa30-Time30(t);
  sumk=sumk+xit(ii)*yy/sigmac0;   
  sumkv=sumkv+xit(ii)^2/sigmac0;
        end
        
        

        end
    end
end

 sumv=(sumkv+eye(1)/(P))\eye(1);


Tr=sumv*sumk;
zeta30=mvnrnd(Tr,sumv); 

%  zeta30=normt_rnd(Tr,sumv,0,9999999);

zeta3s(i,1)=zeta30;

%sample time30
sumtime=0; sumtimev=0;
for t=1:T

    fivt=W19g{t}; omegat=omega0{t};  w5gt=W5g{t}; xit=xi0{t};
    wea1t=wea1{t}; wea2t=wea2{t}; 
  
    for ii=1:n
        for jj=1:n
     
            
        if (ii~=jj)
          r=[fivt(ii,jj),B1(ii,t),B1(jj,t),X31(ii),X31(jj),M(ii,jj),wea1t(ii,jj),wea2t(ii,jj)]; 
        sumtime=sumtime+(w5gt(ii,jj)-r*eta0-omegat(ii)*kappa30-xit(jj)*zeta30)/sigmac0;
        sumtimev=sumtimev+n/sigmac0;
        
        end
        
        

        end
    end
    
     Sigmaf=(0.001+sumtimev)^(-1);
       
        Tf=Sigmaf*sumtime/(sigmac0);
        
        time0=mvnrnd(Tf,Sigmaf);
        
        Time30(t,1)=time0;
    
    
    
end











    




       

end  





 





 results.lambda=lambdas(nomit+1:nit,1);
 results.rho1=rho1s(nomit+1:nit,1);
  results.rho2=rho2s(nomit+1:nit,1);
  results.rho3=rho3s(nomit+1:nit,1);
 results.rho4=rho4s(nomit+1:nit,1);
   results.rho5=rho5s(nomit+1:nit,1);
 results.rho6=rho6s(nomit+1:nit,1);
results.gamma1=gamma1s(nomit+1:nit,1);
results.gamma2=gamma2s(nomit+1:nit,1);
results.gamma3=gamma3s(nomit+1:nit,1);
results.gamma4=gamma4s(nomit+1:nit,1);
  results.beta1=beta1s(nomit+1:nit,1);
  results.beta2=beta2s(nomit+1:nit,1);
  results.beta3=beta3s(nomit+1:nit,1);
results.beta4=beta4s(nomit+1:nit,1);
  results.beta5=beta5s(nomit+1:nit,1);
  results.beta6=beta6s(nomit+1:nit,1);
    results.eta1=eta1s(nomit+1:nit,1);
      results.eta2=eta2s(nomit+1:nit,1);
       results.eta3=eta3s(nomit+1:nit,1);
         results.eta4=eta4s(nomit+1:nit,1);
            results.eta5=eta5s(nomit+1:nit,1);
        results.eta6=eta6s(nomit+1:nit,1);
      results.eta7=eta7s(nomit+1:nit,1);
    results.eta8=eta8s(nomit+1:nit,1);
       results.eta9=eta9s(nomit+1:nit,1);
results.eta10=eta10s(nomit+1:nit,1);

  results.delta1=delta1s(nomit+1:nit,1);
  results.delta2=delta2s(nomit+1:nit,1);
  
  results.delta3=delta3s(nomit+1:nit,1);
  results.delta4=delta4s(nomit+1:nit,1);
  
  results.delta5=delta5s(nomit+1:nit,1);
 results.delta6=delta6s(nomit+1:nit,1);

 

 results.kappa1=kappa1s(nomit+1:nit,1);
 results.kappa2=kappa2s(nomit+1:nit,1);
 results.kappa3=kappa3s(nomit+1:nit,1);
  results.zeta1=zeta1s(nomit+1:nit,1);
 results.zeta2=zeta2s(nomit+1:nit,1);
 results.zeta3=zeta3s(nomit+1:nit,1);
 results.sigmav=sigmavs(nomit+1:nit,1);
 results.sigmae=sigmaes(nomit+1:nit,1);
  results.sigmac=sigmacs(nomit+1:nit,1);

 



end

% =========================================================================
% support functions below
% =========================================================================

function [fvalue ] = likn(Yq, Y10, Y20,Y30,Y40,Y50,Y60, Y11, Y21,Y31,Y41,Y51,Y61, X1q, Xstate,Xtrq, W1q, omegaq,xiq, Time1q,  lambda0,rho0,lambda1,rho1,beta0,kappa10,zeta10,sigmav0)
n=length(Yq);
zzq=[X1q,Xstate,Xtrq];

S0=eye(n)-lambda0*W1q;
S1=eye(n)-lambda1*W1q;

C0=S0*Yq-rho0(1)*Y10-rho0(2)*Y20-rho0(3)*Y30-rho0(4)*Y40-rho0(5)*Y50-rho0(6)*Y60-zzq*beta0-ones(n,1)*Time1q-omegaq*kappa10-xiq*zeta10;
C1=S1*Yq-rho1(1)*Y11-rho1(2)*Y21-rho1(3)*Y31-rho1(4)*Y41-rho1(5)*Y51-rho1(6)*Y61-zzq*beta0-ones(n,1)*Time1q-omegaq*kappa10-xiq*zeta10;

CC0=C0'*C0/(2*sigmav0);
CC1=C1'*C1/(2*sigmav0);

fvalue=(det(S1)/det(S0))*exp(-CC1+CC0);


end
