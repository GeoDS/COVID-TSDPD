clear
rng(101)
data=xlsread('state0108_0415v2.xlsx');

%column 4=state newly confirmed case (NY,OK and MO corrected)
%column 5=state cumulative case 
%column 6=7-day average of new case, (t-1, t-7)
%column 7=7-day average of new case, (t-8, t-14)
%column 8=cross state flow origin 2021
%column 9=cross state flow destination 2021
%column 10=within state flow 2021
%column 11=cross state flow origin 2019
%column 12=cross state flow destination 2019
%column 13=within state flow 2019
%column 14=current temperature
%column 15=1 week lagged temperature
%column 16=2 week lagged temperature
%column 17=3 week lagged temperature
%column 18=4 week lagged temperature
%column 19=wind speed
%column 20= 1 week lagged  wind speed
%column 21= 2 week lagged wind speed
%column 22= 3 week lagged  wind speed
%column 23= 4 week lagged  wind speed
%column 24=precipitation
%column 25=1 week lagged precipitation
%column 26=2 week lagged precipitation
%column 27=3 week lagged precipitation
%column 28=4 week lagged precipitation
%column 29=Binary indicator for zero case day
%column 30=Binary indicator for weekend
%column 31=State population in 2020
%column 32=State population density in 2020
%column 33=5-day average of close contact index, centered at t-5, 2019
%column 34=5-day average of close contact index, centered at t-5, 2021
%column 35=Binary indicator for republican state
%column 36=Median income
%column 37=Percent of male
%column 38=Percent of white
%column 39=Percent of black
%column 40=Percent of asian
%column 41=Percent of Hispanic
%column 42=Percent of population aged 14 and below
%column 43=Percent of population aged 15-64
%column 44=Percent of population aged 65 and above
%column 45=Percent of unemployment
%column 46=Percent of with health  insurance coverage
%column 47=Percentage rank for vaccine
%column 48=percentage rank for population
%column 49=percentage rank for population density
%column 50=percentage rank for age below 14
%column 51=percentage rank for age between 15 and 64
%column 52=percentage rank for age above 64
%column 57=vaccination rate for 1 doses
%column 58=vaccination difference
%column 59=national supply time trend
%column 60=normalized national supply trend
%column 61-158: time trend

yo=data(50:4802,4); % daily new case 
 yc=data(50:4802,5); % daily cumulative case
yl=data(1:4753,4); % lagged daily new case

bo1=data(50:4802,57); %vaccination percentage rate for states
bod1=data(50:4802,58); %change in vaccination rate
bodL=data(1:4753,58); %change in lagged vaccination rate


% average temperature, wind speed, from t-8 to t-14 
 X16=data(50:4802,16);  X21=data(50:4802,21); 

 X15=data(50:4802,15);  X20=data(50:4802,20); 
  % average temperature, wind speed from t-22 to t-28
 X18=data(50:4802,18); X23=data(50:4802,23); 
 
%within state flow
X10=data(50:4802,10);

% dummy for zero case and weekend
X29=data(50:4802,29);
X30=data(50:4802,30);

%vaccination time trend
X60=data(50:4802,60);
%time trend for sdpd
Xtime=data(50:4802,61:158);

%time-invariant state characteristics
X31=data(1:49,31)./1000000; %population in million
X32=data(1:49,32);% population density
X35=data(1:49,35); %republican or not
X36=data(1:49,36); %median income
X37=data(1:49,37); % male ratio
X38=data(1:49,38); % white ratio
X39=data(1:49,39); % black ratio
X40=data(1:49,40); % Asian ratio
X41=data(1:49,41); % Hispanic ratio
X42=data(1:49,42); % pct of below 14;
X43=data(1:49,43); % pct of 15-65
X44=data(1:49,44); % pct of age above 65
X45=data(1:49,45); % pct of unemployment
X46=data(1:49,46); % pct of insurance
X47=data(1:49,47); % pct rank of vaccine
X48=data(1:49,48); % pct rank of population
X49=data(1:49,49); % pct rank of population density



%Xstate=[X33,X34,X37,log(X38),X39/100,X40/100,X41/100,X45/100,X46/100,X47/100,X48/100,X42/100,X43/100];
Xstate=[X31,X32,X35,log(X36),X37/100,X38/100,X39/100,X43/100,X44/100];

n=49; T=97; 

%XB=X49; % percentage quantile of population and population density

Y=zeros(n,T); Yc=zeros(n,T);  YL=zeros(n,T); X1=zeros(n,4,T);  
  Xv=zeros(n,4,T);  Xtr=zeros(n,98,T);


for q=1:T
 yy=yo((q-1)*n+1:q*n,1);

 Y(:,q)=yy;
 
 yyc=yc((q-1)*n+1:q*n,1);

  Yc(:,q)=yyc;
 
 yyl=yl((q-1)*n+1:q*n,1); 

  
 
 YL(:,q)=yyl;
 
 




 
 
% weather controls for  SDPD
%  x16=X16((q-1)*n+1:q*n,1);  X1(:,1,q)=x16;   
% x21=X21((q-1)*n+1:q*n,1);   X1(:,2,q)=x21;    
 x15=X15((q-1)*n+1:q*n,1);  X1(:,1,q)=x15;   
x20=X20((q-1)*n+1:q*n,1);   X1(:,2,q)=x20; 
   x29=X29((q-1)*n+1:q*n,1);    X1(:,3,q)=x29;
   x30=X30((q-1)*n+1:q*n,1);    X1(:,4,q)=x30;    
     
  
   
   % weather controls for vaccination rate
  x18=X18((q-1)*n+1:q*n,1);  Xv(:,1,q)=x18;
 x23=X23((q-1)*n+1:q*n,1);   Xv(:,2,q)=x23; 
   x30=X30((q-1)*n+1:q*n,1);    Xv(:,3,q)=x30; 
   x60=X60((q-1)*n+1:q*n,1);    Xv(:,4,q)=x60; 

% linear and quadratic time trend
Xtr(:,:,q)=Xtime((q-1)*n+1:q*n,:,:);
 
end


M=xlsread('migration_raw.xlsx');
%    Mm=zeros(n,1);
   
 
for i=1:n
    M(i,i)=0;
   
end

  MM=normw(M);

%summ=max(abs(sum(M,2)));

%MM=M./summ;
 
 X32m=median(X32); inpopd=zeros(n,n);
for ii=1:n
    inpopd(ii,ii)=(X32(ii)<X32m);
end

   B1=zeros(n,T); Bd1=zeros(n,T); BdL=zeros(n,T); 
 BBdL=zeros(n,T); 

% B1i=boL(1:49,1);%cumulative vax in period 0
 
 %within state flow
 Fw=zeros(n,T);
 
 for q=1:T
       fw=X10((q-1)*n+1:q*n,1); Fw(:,q)=fw./1000;
    bb1=bo1((q-1)*n+1:q*n,1);  B1(:,q)=bb1;
   bbd1=bod1((q-1)*n+1:q*n,1);  Bd1(:,q)=bbd1;
   bbdL=bodL((q-1)*n+1:q*n,1);    BdL(:,q)=bbdL;
  BBdL(:,q)=MM*bbdL;
%    Br(:,:,q)=[bb1, XB,ones(n,1)];

 end
   
   

  load StateflowV_weights.mat;
 

%date 0 is 01-08, data 1 is 01-09;
W5=cell(T,1); W5g=cell(T,1);


ww1=xlsread('2021-01-08.xlsx'); ww2=xlsread('2021-01-09.xlsx'); ww3=xlsread('2021-01-10.xlsx');
ww4=xlsread('2021-01-11.xlsx'); ww5=xlsread('2021-01-12.xlsx'); ww6=xlsread('2021-01-13.xlsx');
ww7=xlsread('2021-01-14.xlsx'); ww8=xlsread('2021-01-15.xlsx'); ww9=xlsread('2021-01-16.xlsx');
ww10=xlsread('2021-01-17.xlsx'); ww11=xlsread('2021-01-18.xlsx'); ww12=xlsread('2021-01-19.xlsx');
ww13=xlsread('2021-01-20.xlsx'); ww14=xlsread('2021-01-21.xlsx'); ww15=xlsread('2021-01-22.xlsx');
ww16=xlsread('2021-01-23.xlsx'); ww17=xlsread('2021-01-24.xlsx'); ww18=xlsread('2021-01-25.xlsx');
ww19=xlsread('2021-01-26.xlsx'); ww20=xlsread('2021-01-27.xlsx'); ww21=xlsread('2021-01-28.xlsx');
ww22=xlsread('2021-01-29.xlsx'); ww23=xlsread('2021-01-30.xlsx'); ww24=xlsread('2021-01-31.xlsx');
ww25=xlsread('2021-02-01.xlsx'); ww26=xlsread('2021-02-02.xlsx'); ww27=xlsread('2021-02-03.xlsx');
ww28=xlsread('2021-02-04.xlsx');

W5{1}=ww2./1000; W5{2}=ww3./1000; W5{3}=ww4./1000; W5{4}=ww5./1000; W5{5}=ww6./1000; W5{6}=ww7./1000;
W5{7}=ww8./1000; W5{8}=ww9./1000; W5{9}=ww10./1000; W5{10}=ww11./1000; W5{11}=ww12./1000; W5{12}=ww13./1000; 
W5{13}=ww14./1000; W5{14}=ww15./1000; W5{15}=ww16./1000; W5{16}=ww17./1000; W5{17}=ww18./1000; W5{18}=ww19./1000;
 W5{19}=ww20./1000; W5{20}=ww21./1000; W5{21}=ww22./1000; W5{22}=ww23./1000; W5{23}=ww24./1000; W5{24}=ww25./1000;
 W5{25}=ww26./1000; W5{26}=ww27./1000; W5{27}=ww28./1000;

 W5g{1}=ww2./1000; W5g{2}=ww3./1000; W5g{3}=ww4./1000; W5g{4}=ww5./1000; W5g{5}=ww6./1000; W5g{6}=ww7./1000;
W5g{7}=ww8./1000; W5g{8}=ww9./1000; W5g{9}=ww10./1000; W5g{10}=ww11./1000; W5g{11}=ww12./1000; W5g{12}=ww13./1000; 
W5g{13}=ww14./1000; W5g{14}=ww15./1000; W5g{15}=ww16./1000; W5g{16}=ww17./1000; W5g{17}=ww18./1000; W5g{18}=ww19./1000;
 W5g{19}=ww20./1000; W5g{20}=ww21./1000; W5g{21}=ww22./1000; W5g{22}=ww23./1000; W5g{23}=ww24./1000; W5g{24}=ww25./1000;
 W5g{25}=ww26./1000; W5g{26}=ww27./1000; W5g{27}=ww28./1000;


 for t=1:T-27
    
 ww=WWv{t};

 W5{t+27}=ww./1000;
 W5g{t+27}=ww./1000;
 
end



rsum=zeros(T,1); 

for t=1:T
    rsum(t,1)=max(sum(W5{t},2));
   
end

factor=max(max(rsum));

%scalar normalization for spatial weights on time t in SDPD
for t=1:T
   
 ww=W5{t};
ww=ww./factor;
   W5{t}=ww;

   
end


% 
rsum=zeros(T,1);  
for t=1:T
 rsum(t,1)=max(sum(W5{t},2));   
 end

rs=max(rsum);

% 2019 un-row normalized spatial weights
load Stateflow2019_weights.mat;
W19g=cell(T,1);


 ww2=xlsread('2019-01-09.xlsx'); ww3=xlsread('2019-01-10.xlsx');
ww4=xlsread('2019-01-11.xlsx'); ww5=xlsread('2019-01-12.xlsx'); ww6=xlsread('2019-01-13.xlsx');
ww7=xlsread('2019-01-14.xlsx'); ww8=xlsread('2019-01-15.xlsx'); ww9=xlsread('2019-01-16.xlsx');
ww10=xlsread('2019-01-17.xlsx'); ww11=xlsread('2019-01-18.xlsx'); ww12=xlsread('2019-01-19.xlsx');
ww13=xlsread('2019-01-20.xlsx'); ww14=xlsread('2019-01-21.xlsx'); ww15=xlsread('2019-01-22.xlsx');
ww16=xlsread('2019-01-23.xlsx'); ww17=xlsread('2019-01-24.xlsx'); ww18=xlsread('2019-01-25.xlsx');
ww19=xlsread('2019-01-26.xlsx'); ww20=xlsread('2019-01-27.xlsx'); ww21=xlsread('2019-01-28.xlsx');
ww22=xlsread('2019-01-29.xlsx'); ww23=xlsread('2019-01-30.xlsx'); ww24=xlsread('2019-01-31.xlsx');
ww25=xlsread('2019-02-01.xlsx'); ww26=xlsread('2019-02-02.xlsx'); ww27=xlsread('2019-02-03.xlsx');
ww28=xlsread('2019-02-04.xlsx');

W19g{1}=ww2./1000; W19g{2}=ww3./1000; W19g{3}=ww4./1000; W19g{4}=ww5./1000; W19g{5}=ww6./1000; W19g{6}=ww7./1000;
W19g{7}=ww8./1000; W19g{8}=ww9./1000; W19g{9}=ww10./1000; W19g{10}=ww11./1000; W19g{11}=ww12./1000; W19g{12}=ww13./1000; 
W19g{13}=ww14./1000; W19g{14}=ww15./1000; W19g{15}=ww16./1000; W19g{16}=ww17./1000; W19g{17}=ww18./1000; W19g{18}=ww19./1000;
 W19g{19}=ww20./1000; W19g{20}=ww21./1000; W19g{21}=ww22./1000; W19g{22}=ww23./1000; W19g{23}=ww24./1000; W19g{24}=ww25./1000;
 W19g{25}=ww26./1000; W19g{26}=ww27./1000; W19g{27}=ww28./1000;





for t=1:T-27
   w19=WWp{t}; 
   W19g{t+27}=w19./1000; 
end


 %weather difference control for gravity
 wea1=cell(T,1);  wea2=cell(T,1); wea3=cell(T,1);
for t=1:T
    XX1=zeros(n,n); XX2=zeros(n,n); XX3=zeros(n,n);
    XX4=zeros(n,n); XX5=zeros(n,n);
    for i=1:n
        for j=1:n
            if (i~=j)
               XX1(i,j)=abs(X1(i,1,t)-X1(j,1,t));
               XX2(i,j)=abs(X1(i,2,t)-X1(j,2,t));
               XX3(i,j)=abs(X1(i,3,t)-X1(j,3,t));
             
            end
        end
    end
    
    wea1{t}=XX1; wea2{t}=XX2; wea3{t}=XX3;

  
    
end


% store national newly confirmed and cumulative cases
Yn=zeros(T,1); % national new cases
Vn=zeros(T,1); %national facutal cumulative vax rate
Vnewn=zeros(T,1); %national factual new vax rate
Ync=zeros(T,1);  %national cumulative cases
Ync(1,1)=sum(Y(:,1))+sum(YL(:,1));
Ycstate=zeros(n,T); %state cumulative case
Ycstate(:,1)=Y(:,1)+YL(:,1);
popweight=X31./(sum(X31));

for t=1:T
 Yn(t,1)=sum(Y(:,t)); 
vcp=0;  vnewcp=0;
 for ii=1:n
   vcp=vcp+popweight(ii)*B1(ii,t); 
   vnewcp=vnewcp+popweight(ii)*Bd1(ii,t); 
 end
Vn(t,1)=vcp;
Vnewn(t,1)=vnewcp;
 end



for t=2:T
    Ync(t,1)=Ync(t-1,1)+Yn(t,1);
    Ycstate(:,t)=Ycstate(:,t-1)+Y(:,t);
end

B11=sort(B1);
lgamma=min(B11(round(length(B11(:,2))*0.05),:));
ugamma=max(B11(round(length(B11(:,2))*0.95),:));


nit=25000; br=0.2;
nomit=nit*br;

indiv10=zeros(n,1); indiv20=zeros(n,1); 
Time10=zeros(T,1); 
for q=3:T
Time10(q,1)=0.05;    
end
sigmav0=1; sigmae0=1; sigmac0=1; 
Time20=Time10;  Time30=Time10;

lambda0=0.2;

rho0=[0.5;0.2;0;0.5;0.3;0.1];
gamma0=[3;10; 6; 11]; 
% gamma0=[3;5];

kappa0=0.5*ones(3,1); zeta0=0.5*ones(3,1); 
delta0=0.5*ones(15,1);
beta0=0.5*ones(111,1); 
eta0=0.5*ones(8,1);


sce=1; % 1 for muting vaccination effect on within state transmission
      % 2 for muting vaccination effect on cross state flow
       % 3 for muting all vaccination channels






results=mcmctsdpdbaselinei(nit,br, n,T,Y,YL,X1,B1,Bd1,BdL,BBdL,Xv,X31,Xstate,Xtr,M,W5,inpopd,W5g,W19g,wea1,wea2,Time10,Time30,lambda0, rho0,gamma0,lgamma,ugamma,beta0, delta0,eta0,kappa0,zeta0,sigmav0,sigmae0,sigmac0,rs);
lambdas=results.lambda;
rho1s=results.rho1;
rho2s=results.rho2;
rho3s=results.rho3;
rho4s=results.rho4;
rho5s=results.rho5;
rho6s=results.rho6;
gamma1s=results.gamma1;
gamma2s=results.gamma2;
gamma3s=results.gamma3;
gamma4s=results.gamma4;
beta1s=results.beta1;  % temperature
beta2s=results.beta2;  % wind speed
beta3s=results.beta3;
beta4s=results.beta4;
beta5s=results.beta5;
beta6s=results.beta6;
delta1s=results.delta1;
delta2s=results.delta2;
delta3s=results.delta3;
delta4s=results.delta4;
delta5s=results.delta5;
delta6s=results.delta6;

kappa1s=results.kappa1;
kappa2s=results.kappa2;
kappa3s=results.kappa3;
zeta1s=results.zeta1;
zeta2s=results.zeta2;
zeta3s=results.zeta3;


eta1s=results.eta1; % 2019 cross-state flow
eta2s=results.eta2;  % destination vaccination influence
eta3s=results.eta3; % origin vaccination influence
eta4s=results.eta4;
eta5s=results.eta5;
eta6s=results.eta6;
eta7s=results.eta7; %temperature difference
eta8s=results.eta8; % wind difference



  


%   histogram(lambdas,'Normalization', 'pdf');
% hold on;
% ksdensity(lambdas);
% hold on;
% xline(0.2, 'r--', 'LineWidth', 2, 'Label', 'True');
% 
%   histogram(zeta1s,'Normalization', 'pdf');
% hold on;
% ksdensity(zeta1s);
% hold on;
% xline(mean(zeta1s), 'r--', 'LineWidth', 2, 'Label', 'True');
% 
% 
% 
% 
% histogram(rho1s,'Normalization', 'pdf');
% hold on;
% ksdensity(rho1s);
% xline(mean(rho1s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% 
% histogram(rho2s,'Normalization', 'pdf');
% hold on;
% ksdensity(rho2s);
% xline(mean(rho2s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(rho3s,'Normalization', 'pdf');
% hold on;
% ksdensity(rho3s);
% xline(mean(rho3s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(rho4s,'Normalization', 'pdf');
% hold on;
% ksdensity(rho4s);
% xline(mean(rho4s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% 
% histogram(rho5s,'Normalization', 'pdf');
% hold on;
% ksdensity(rho5s);
% xline(mean(rho5s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(rho6s,'Normalization', 'pdf');
% hold on;
% ksdensity(rho6s);
% xline(mean(rho6s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(gamma1s,'Normalization', 'pdf');
% hold on;
% ksdensity(gamma1s);
% xline(mean(gamma1s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(gamma2s,'Normalization', 'pdf');
% hold on;
% ksdensity(gamma2s);
% xline(mean(gamma2s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(gamma3s,'Normalization', 'pdf');
% hold on;
% ksdensity(gamma3s);
% xline(mean(gamma3s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% 
% histogram(gamma4s,'Normalization', 'pdf');
% hold on;
% ksdensity(gamma4s);
% xline(mean(gamma4s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% 
% 
% histogram(kappa1s,'Normalization', 'pdf');
% hold on;
% ksdensity(kappa1s);
% xline(mean(kappa1s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(kappa2s,'Normalization', 'pdf');
% hold on;
% ksdensity(kappa2s);
% xline(mean(kappa2s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(kappa3s,'Normalization', 'pdf');
% hold on;
% ksdensity(kappa3s);
% xline(mean(kappa3s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% 
% 
% histogram(zeta1s,'Normalization', 'pdf');
% hold on;
% ksdensity(zeta1s);
% xline(mean(zeta1s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(zeta2s,'Normalization', 'pdf');
% hold on;
% ksdensity(zeta2s);
% xline(mean(zeta2s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% histogram(zeta3s,'Normalization', 'pdf');
% hold on;
% ksdensity(zeta3s);
% xline(mean(zeta3s), 'r--', 'LineWidth', 2, 'Label', 'posterior mean');
% 
% 
% 
% 
% 
% 
% 
% %  % store by state CF new case and cumulative case
% % filename='CF_case_s3.xlsx';
% % xlswrite(filename, Ymean);
% % 
% %  filename='CF_case_ub_s3.xlsx';
% % xlswrite(filename, Yub);
% % 
% %  filename='CF_case_lb_s3.xlsx';
% % xlswrite(filename, Ylb);
% % %
% % filename='Factual_case.xlsx';
% % xlswrite(filename, Y);
% % % 
% % % 
% %  filename='CF_ccase_s3.xlsx';
% % xlswrite(filename, Ycmean);
% % 
% %  filename='CF_ccase_ub_s3.xlsx';
% % xlswrite(filename, Ycub);
% % 
% %  filename='CF_ccase_lb_s3.xlsx';
% % xlswrite(filename, Yclb);
% % 
% % filename='Factual_ccase.xlsx';
% % Ycend=Ycstate(:,T);
% % xlswrite(filename, Ycend);
% % 
% % 
% % 
% % %********************************
% % % plot national trend of CF new case
% % % ******************************
% %   tn = datetime(2021,1,8) + caldays(0:T-1);
% %  dn=datenum(tn);
% %   hperc=0.95;
% % 
% %  Ycm=mean(Ycs,2);  
% % 
% %  YcsU=zeros(T,1);  YcsL=zeros(T,1);
% % 
% % 
% %  for t=1:T
% %      bound=hpdi(Ycs(t,:)',hperc);
% %      YcsU(t,1)=bound(1);  YcsL(t,1)=bound(2);
% %  end   
% 
% % 
% % 
% % 
% %    plot(dn,Yn(1:T),'k');
% %  datetick('x','yyyy-mm-dd');
% %  hold on;
% %  plot(dn,Ycm(1:T),'--k');
% %  datetick('x','yyyy-mm-dd');  
% %  hold on;
% %  plot(dn,YcsU,'o');
% %  hold on;
% %  plot(dn,YcsL,'*');
% %  
%   filename='CF_nationalnewcase_s3.xlsx';
% xlswrite(filename, Ycm);
% 
%  filename='CF_nationalnewcase_ub_s3.xlsx';
% xlswrite(filename, YcsU);
% 
%  filename='CF_nationalnewcase_lb_s3.xlsx';
% xlswrite(filename, YcsL);
% % %
% filename='Factual_nationalnewcase.xlsx';
% xlswrite(filename, Yn);
% 
% 
% 
% %*****************************
% % plot national trend of CF cumulative case
%  %********************************** 
% 
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
% % 
% % 
% % 
% %    plot(dn,Ycumn(1:T),'k');
% %  datetick('x','yyyy-mm-dd');
% %  hold on;
% %  plot(dn,Yccumn(1:T),'--k');
% %  datetick('x','yyyy-mm-dd');  
% %  hold on;
% %  plot(dn,YctU,'o');
% %  hold on;
% %  plot(dn,YctL,'*');
% %  
% %  
%      filename='CF_nationalcumulativecase_s3.xlsx';
% xlswrite(filename, Yccumn);
% 
%  filename='CF_nationalcumulativecase_ub_s3.xlsx';
% xlswrite(filename, YctU);
% 
%  filename='CF_nationalcumulativecase_lb_s3.xlsx';
% xlswrite(filename, YctL);
% % %
% filename='Factual_nationalcumulativecase.xlsx';
% xlswrite(filename, Ync);
% 
% 
% 
% 
% 
% 
% 
