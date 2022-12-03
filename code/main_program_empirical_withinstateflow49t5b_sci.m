clear
rng(101)
data=xlsread('state0205_0415v13.xlsx');

%column 4=state newly confirmed case (NY,OK and MO corrected)
%column 5=state cumulative case 
%column 6=state newly confirmed case  7 day average
%column 7=state cumulative 7 day average
%column 8=state newly confirmed case per capita
%column 9=state cumulative case per capita
%column 10=state new case 7 day average per capita
%column 11=state cumulative case 7 day average per capita
%column 12=weekly growth rate by victor
%column 13=vaccination rate centered at t-19
%column 14=IV for vaccination
%column 15=cross state flow 2021
%column 16=within state flow 2021
%column 17=cross state flow 2019
%column 18=within state flow 2019
%column 19=state within/cross county ratio 2021
%column 20=state within/cross county ratio 2019
%column 21=current temperature
%column 22=1 week lagged temperature
%column 23=2 week lagged temperature
%column 24=3 week lagged temperature
%column 25=4 week lagged temperature
%column 26=current wind speed
%column 27=1 week lagged wind speed
%column 28=2 week lagged wind speed
%column 29=3 week lagged wind speed
%column 30=4 week lagged wind speed
%column 31=maximum wind speed
%column 32= 1 week lagged maximum wind speed
%column 33= 2 week lagged maximum wind speed
%column 34= 3 week lagged maximum wind speed
%column 35= 4 week lagged maximum wind speed
%column 36=precipitation
%column 37=1 week lagged precipitation
%column 38=2 week lagged precipitation
%column 39=3 week lagged precipitation
%column 40=4 week lagged precipitation
%column 46=non reporting dummy
%column 47=state population
%column 48=2020 close contact index
%column 49=2021 close contact index
%column 51=updated vax completion rate from CDC
%column 52=corrected vax completion rate from github


yo=data(50:3430,4); % daily new case 7 day average
yl=data(1:3381,4); % lagged daily confirm case 7 day average
yc=data(50:3430,5); %state cumulative case 7 day average


bo1=data(50:3430,52); %vaccination percentage rate for states
boL=data(1:3381,52); % lagged vaccination rate

% average temperature, wind speed, maximum wind speed, and precipitation
% from t-8 to t-14 
 X23=data(50:3430,23);  X28=data(50:3430,28); X33=data(50:3430,33); X38=data(50:3430,38);

 % average temperature, wind speed, maximum wind speed, and precipitation
% from t-22 to t-28
X25=data(50:3430,25); X30=data(50:3430,30); X35=data(50:3430,35); X40=data(50:3430,40);

% IV for vaccination
X14=data(50:3430,14);

%2021 within/cross state ratio
X16=data(50:3430,16);

%2019 within/cross state ratio
X18=data(50:3430,18);

%non reporting dummy
X46=data(50:3430,46);

%population in million
X47=data(1:49,47)./1000000;
popweight=X47./sum(X47);



n=49; T=69; 

X47m=median(X47);
inpop=zeros(n,n);
for ii=1:n
    inpop(ii,ii)=(X47(ii)<X47m);
end


Y=zeros(n,T);   YL=zeros(n,T); X1=zeros(n,4,T);   nr=zeros(n,1,T);
 ew=zeros(n,4,T); Yc=zeros(n,T);

for q=1:T
 yy=yo((q-1)*n+1:q*n,1);
 
Y(:,q)=yy;
 
 yyc=yc((q-1)*n+1:q*n,1);

 Yc(:,q)=yyc;
 
 yyl=yl((q-1)*n+1:q*n,1); 

 YL(:,q)=yyl;
 
 
% weather controls for within/cross county ratio and SDPD
      x23=X23((q-1)*n+1:q*n,1);  X1(:,1,q)=x23;
 x28=X28((q-1)*n+1:q*n,1);    X1(:,2,q)=x28; 
   x33=X33((q-1)*n+1:q*n,1);   X1(:,3,q)=x33; 
   x38=X38((q-1)*n+1:q*n,1);    X1(:,4,q)=x38;
   x46=X46((q-1)*n+1:q*n,1);    nr(:,1,q)=x46;
   
 % weather controls for vaccination rate
   x25=X25((q-1)*n+1:q*n,1);  ew(:,1,q)=x25;
  x30=X30((q-1)*n+1:q*n,1);   ew(:,2,q)=x30;
   x35=X35((q-1)*n+1:q*n,1);   ew(:,3,q)=x35; 
     x40=X40((q-1)*n+1:q*n,1);    ew(:,4,q)=x40;     
end

B1=zeros(n,T);  E1=zeros(n,T); Fw=zeros(n,T);
Fw19=zeros(n,T); BL=zeros(n,T); Ci=zeros(n,T);
Bp=zeros(n,T);   

   for q=1:T
    bb1=bo1((q-1)*n+1:q*n,1);  B1(:,q)=bb1;
   bbL=boL((q-1)*n+1:q*n,1);    BL(:,q)=bbL;
   ee1=X14((q-1)*n+1:q*n,1);  E1(:,q)=ee1;
    fw=X16((q-1)*n+1:q*n,1);  Fw(:,q)=fw./1000;
 fw19=X18((q-1)*n+1:q*n,1);   Fw19(:,q)=fw19./1000;

   end
   
M=xlsread('fb.xlsx');
 
for i=1:n
  M(i,i)=0;  
end

MM=normw(M);


  load StateflowV_weights.mat;
 
W5=cell(T,1); W5L=cell(T,1);  
W5g=cell(T,1);

 for t=1:T
    
 ww=WWv{t+1};

 W5{t}=ww./1000;
 W5g{t}=ww./1000;
 
end

for t=1:T
  w1=WWv{t};
 
 W5L{t}=w1./1000;
end

rsum=zeros(T,1);  rsuml=zeros(T,1);

for t=1:T
    rsum(t,1)=max(sum(W5{t},2));
    rsuml(t,1)=max(sum(W5L{t},2));
end

factor=max(max(rsum),max(rsuml));

for t=1:T
    %row-normalized population flow from Wuhan
 ww=W5{t};
ww=ww./factor;
   W5{t}=ww;

   
 end

for t=1:T

   w1=W5L{t};

w1=w1./factor;
      
    W5L{t}=w1;
end
% 
rsum=zeros(T,1);  rsuml=zeros(T,1);
for t=1:T
 rsum(t,1)=max(sum(W5{t},2));   
 rsuml(t,1)=max(sum(W5L{t},2)); 
end

rs=max(max(rsum),max(rsuml));


load Stateflow2019_weights.mat;
W19g=cell(T,1);
for t=1:T
   w19=WWp{t+1}; 
   W19g{t}=w19./1000; 
end





%weather difference control for gravity
wea1=cell(T,1);  wea2=cell(T,1); wea3=cell(T,1); wea4=cell(T,1);
for t=1:T
    XX1=zeros(n,n); XX2=zeros(n,n); XX3=zeros(n,n); XX4=zeros(n,n);
  
    for i=1:n-1
        for j=1:n-1
            if (i~=j)
               XX1(i,j)=abs(X1(i,1,t)-X1(j,1,t));
               XX2(i,j)=abs(X1(i,2,t)-X1(j,2,t));
               XX3(i,j)=abs(X1(i,3,t)-X1(j,3,t));
               XX4(i,j)=abs(X1(i,4,t)-X1(j,4,t));
             
            end
        end
    end
    
    wea1{t}=XX1; wea2{t}=XX2; wea3{t}=XX3; wea4{t}=XX4;
  
    
end


% store national newly confirmed and cumulative cases
Yn=zeros(T,1);  Vn=zeros(T,1);
Ynstate=zeros(n,T); %state cumulative case
Ynstate(:,1)=Y(:,1)+YL(:,1);

for t=1:T
 Yn(t,1)=sum(Y(:,t));  
 vcp=0;
 for ii=1:n
   vcp=vcp+popweight(ii)*B1(ii,t)/n; 
 end

Vn(t,1)=vcp;
end

for t=2:T
 Ynstate(:,t)=Ynstate(:,t-1)+Y(:,t);   
end





B11=sort(B1);
lgamma=min(B11(round(length(B11(:,2))*0.1),:));
ugamma=max(B11(round(length(B11(:,2))*0.9),:));


nit=25000; br=0.5;
nomit=nit*br;





indiv10=zeros(n,1); indiv20=zeros(n,1); indiv30=zeros(n,1);
Time10=zeros(T,1);
for q=3:T
Time10(q,1)=0.05;    
end
sigmav0=1; sigmae0=1; sigmat0=1; 
Time20=Time10; Time30=Time20;

beta0=0.5*ones(6,1); 
% beta0=0.5*ones(7,1); 
kappa10=0.5; kappa20=0.5; kappa30=0.5;
kappa40=0.5; kappa50=0.5;


sce=5; % 2 for muting vaccination effect on within state transmission
       % 3 for muting vaccination effect on within state flow
       % 4 for muting vaccination effect on cross state flow
       % 5 for muting all vaccination channels

   
sigmam0=1;
lambda10=0.2;  lambda0=lambda10;%

rho10=0.2; rho20=0;
rho30=0.6; rho40=0.4; %  gamma0=2; rho40=0;
gamma0=[0.6;2.8];
 rho0=[rho10;rho20;rho30;rho40]; 
mu0=0; eta0=ones(10,1);
delta0=0.3*ones(6,1);
phi0=0.5*ones(6,1);
kappa0=0.5*ones(4,1); zeta0=0.5*ones(4,1); 



results=mcmcwfbd8b(nit,br, n,T,Y,YL,nr, X1,ew,B1,BL,E1, Fw,Fw19,M,MM,W5,W5L,W5g,W19g,inpop,wea1,wea2,wea3,wea4,lambda0, rho0,mu0,gamma0,lgamma,ugamma,beta0, delta0,phi0,eta0,kappa0,zeta0,sigmav0,sigmae0,sigmat0,sigmam0,X47,indiv10,indiv20,indiv30,Time10,Time30,rs,sce);

lambdas=results.lambda;
   rho1s=results.rho1;
    rho2s=results.rho2;
       rho3s=results.rho3;
    rho4s=results.rho4;
    mus=results.mu;
    gamma1s=results.gamma1;
     gamma2s=results.gamma2;
    beta1s=results.beta1;
   beta2s=results.beta2;
   beta3s=results.beta3;
   beta4s=results.beta4;
   beta5s=results.beta5;
   beta6s=results.beta6;
    delta1s=results.delta1;
    delta2s=results.delta2;
    delta3s=results.delta3;
    delta4s=results.delta4;
    delta5s=results.delta5;
    phi1s=results.phi1;
      phi2s=results.phi2; 
      phi3s=results.phi3;
      phi4s=results.phi4;
      phi5s=results.phi5;
      sigmavs=results.sigmav;
      sigmams=results.sigmam;
      sigmats=results.sigmat;
 

   eta1s=results.eta1;
   eta2s=results.eta2;
   eta3s=results.eta3;
   eta4s=results.eta4;
   eta5s=results.eta5;
     eta6s=results.eta6;
   eta7s=results.eta7;
   
  
%   Ycs=results.Yc;
% Vcs=results.Vc;
%   Ymean=results.Ymean;
%   Yub=results.Yub;
%   Ylb=results.Ylb;
%   
%     Ycmean=results.Ycmean;
%   Ycub=results.Ycub;
%   Yclb=results.Yclb;
%   
%     Vmean=results.Vmean;
%   Vub=results.Vub;
%   Vlb=results.Vlb;
%   
%   
%   Wp=results.Wp;
%   outflowm=results.outflowm;
%   outflowu=results.outflowu;
% outflowl=results.outflowl;
% 
%   inflowm=results.inflowm;
%   inflowu=results.inflowu;
% inflowl=results.inflowl;

% 

% % 
%  filename='CF_network_flow12.xlsx';
% xlswrite(filename, Wp);
%  filename='Factual_network.xlsx';
% xlswrite(filename, Wf);


%  filename='CF_case_s5.xlsx';
% xlswrite(filename, Ymean);
% 
%  filename='CF_case_ub_s5.xlsx';
% xlswrite(filename, Yub);
% 
%  filename='CF_case_lb_s5.xlsx';
% xlswrite(filename, Ylb);
% %
% filename='Factual_case.xlsx';
% xlswrite(filename, Y);
% 
% 
%  filename='CF_ccase_s5.xlsx';
% xlswrite(filename, Ycmean);
% 
%  filename='CF_ccase_ub_s5.xlsx';
% xlswrite(filename, Ycub);
% 
%  filename='CF_ccase_lb_s5.xlsx';
% xlswrite(filename, Yclb);
% filename='Factual_ccase.xlsx';
% Ycend=Yc(:,T);
% xlswrite(filename, Ycend);
% 



% % 
%  filename='CF_outflow_s5.xlsx';
% xlswrite(filename, outflowm);
% 
%  filename='CF_outflow_ub_s5.xlsx';
% xlswrite(filename, outflowu);
% 
%  filename='CF_outflow_lb_s5.xlsx';
% xlswrite(filename, outflowl);
% filename='Factual_outflow.xlsx';
% xlswrite(filename, outflow);
% 
%  filename='CF_inflow_s5.xlsx';
% xlswrite(filename, inflowm);
% 
%  filename='CF_inflow_ub_s5.xlsx';
% xlswrite(filename, inflowu);
% 
%  filename='CF_inflow_lb_s5.xlsx';
% xlswrite(filename, inflowl);
% filename='Factual_inflow.xlsx';
% xlswrite(filename, inflow);
% 



% % 
% 
% % % 
%  filename='CF_vax_flow12.xlsx';
% xlswrite(filename, Vmean);
% 
%  filename='CF_vax_ub_flow12.xlsx';
% xlswrite(filename, Vub);
% 
%  filename='CF_vax_lb_flow12.xlsx';
% xlswrite(filename, Vlb);
% 
%   filename='Factual_vax.xlsx';
% xlswrite(filename, B1);

  
 

% hperc=0.95;
% % plot newly confirmed cases
% tn = datetime(2021,2,6) + caldays(0:T-1);
% dn=datenum(tn);
% 
% %**********************************************************
%  % generate counterfactual for newly confirm cases 
%  %*********************************************************
%  
% Vcm=mean(Vcs,2);  
% VcsU=zeros(T,1);  VcsL=zeros(T,1);
%  
%  
%  for t=1:T
%      bound=hpdi(Vcs(t,:)',hperc);
%      VcsU(t,1)=bound(1);  VcsL(t,1)=bound(2);
%  end   
%  
% 
% 
% 
%    plot(dn,Vn(1:T),'k');
%  datetick('x','yyyy-mm-dd');
%  hold on;
%  plot(dn,Vcm(1:T),'--k');
%  datetick('x','yyyy-mm-dd');  
%  hold on;
%  plot(dn,VcsU,'o');
%  hold on;
%  plot(dn,VcsL,'*');

  
%   hperc=0.95;
% % plot newly confirmed cases
% tn = datetime(2021,2,6) + caldays(0:T-1);
% dn=datenum(tn);
% 
% %**********************************************************
%  % generate counterfactual for newly confirm cases 
%  %*********************************************************
%  
% 
% 
% 
% Ycm=mean(Ycs,2);  
% 
%  
%  
%  
%  
%  YcsU=zeros(T,1);  YcsL=zeros(T,1);
%  
%  
%  for t=1:T
%      bound=hpdi(Ycs(t,:)',hperc);
%      YcsU(t,1)=bound(1);  YcsL(t,1)=bound(2);
%  end   
%  
% 
% 
% 
%    plot(dn,Yn(1:T),'k');
%  datetick('x','yyyy-mm-dd');
%  hold on;
%  plot(dn,Ycm(1:T),'--k');
%  datetick('x','yyyy-mm-dd');  
%  hold on;
%  plot(dn,YcsU,'o');
%  hold on;
%  plot(dn,YcsL,'*');
%  
%   filename='CF_nationalnewcase_s5.xlsx';
% xlswrite(filename, Ycm);
% 
%  filename='CF_nationalnewcase_ub_s5.xlsx';
% xlswrite(filename, YcsU);
% 
%  filename='CF_nationalnewcase_lb_s5.xlsx';
% xlswrite(filename, YcsL);
% % %
% filename='Factual_nationalnewcase.xlsx';
% xlswrite(filename, Yn);
% %  
% %  
% % 
% %    
% % %**********************************************************
% %  % generate counterfactual for cumulative cases
% %  %*********************************************************
% %  
% %  % compute the actual national cumulative case 
% Ycumn=zeros(T,1); 
% for t=1:T
%   Ycumn(t,1)=sum(Yc(:,t));  
% end
% 
% %derive the mcmc draws of the predicted counterfactual national cumulative case
% Yct=zeros(T,nit-nomit);
%  Yct(1,:)=Ycs(1,:)+ones(1,nit-nomit)*sum(YL(:,1));
%  for t=2:T
%    Yct(t,:)=Yct(t-1,:)+Ycs(t,:);
%     
%  end
% 
%  %compute the mean counterfactual national cumulative case, and CI
% Yccumn=zeros(T,1);
% 
% Ycm=mean(Ycs,2);  
% 
% Yccumn(1,1)=Ycm(1,1)+sum(YL(:,1));
% for t=2:T
%    Yccumn(t,1)=Yccumn(t-1,1)+Ycm(t,1);
%     
% end
% 
%  
%  
%  
%  
%  YctU=zeros(T,1);  YctL=zeros(T,1);
%  
%  
%  for t=1:T
%      bound=hpdi(Yct(t,:)',hperc);
%      YctU(t,1)=bound(1);  YctL(t,1)=bound(2);
%  end   
%  
% 
% 
% 
%    plot(dn,Ycumn(1:T),'k');
%  datetick('x','yyyy-mm-dd');
%  hold on;
%  plot(dn,Yccumn(1:T),'--k');
%  datetick('x','yyyy-mm-dd');  
%  hold on;
%  plot(dn,YctU,'o');
%  hold on;
%  plot(dn,YctL,'*');
%  
%  
%      filename='CF_nationalcumulativecase_s5.xlsx';
% xlswrite(filename, Yccumn);
% 
%  filename='CF_nationalcumulativecase_ub_s5.xlsx';
% xlswrite(filename, YctU);
% 
%  filename='CF_nationalcumulativecase_lb_s5.xlsx';
% xlswrite(filename, YctL);
% % %
% filename='Factual_nationalcumulativecase.xlsx';
% xlswrite(filename, Ycumn);