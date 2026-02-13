close all;
clear all;
K=150; 

%{
figure
hold on
term_10 = importdata('term_tracking_10_long.tab');
endpt_25 = importdata('term_tracking_25_endpt_long.tab');
term_15 = importdata('term_tracking_15_long.tab');
plot(0:K-1,endpt_25(1:K,2),0:K-1,term_15(1:K,2),'-.',0:K-1,term_10(1:K,2),'--')
stairs([0,50,100,150],[.41,.36,.38,.38],'Linewidth',1.75)
legend('N=25, endpt','N=15, term con','N=10, term con','Terminal Region')
xlabel('k')
ylabel('$|x_k - x_{ss}|$')


figure
hold on
term_10_unc = importdata('term_tracking_10_unc_long.tab');
endpt_25_unc = importdata('term_tracking_25_endpt_unc_long.tab');
term_15_unc = importdata('term_tracking_15_unc_long.tab');
plot(0:K-1,endpt_25_unc(1:K,2),0:K-1,term_15_unc(1:K,2),'-.',0:K-1,term_10_unc(1:K,2),'--')
stairs([0,50,100,150],[.41,.36,.38,.38],'Linewidth',1.75)
legend('N=25, endpt','N=15, term con','N=10, term con','Terminal Region')
xlabel('k')
ylabel('$|x_k - x_{ss}|$')

nomtime=[ mean(endpt_25(1:K,4)), mean(term_15(1:K,4)) , mean(term_10(1:K,4)) ] 
unctime=[ mean(endpt_25_unc(1:K,4)), mean(term_15_unc(1:K,4)) , mean(term_10_unc(1:K,4)) ]
%}

lecss = [-0.222801341558419*ones(K/3,1)  ;   -0.511719641758650*ones(K/3,1) ; -0.520230985947510*ones(K/3,1)] ; 

ec_nu = importdata('ec_nu.tab');
tracking_nu = importdata('tracking_nu.tab');
reg_nu = importdata('reg_nu.tab');
scp99_nu = importdata('scp99_nu.tab');
sc2_nu = importdata('sc2p5_nu.tab');
figure
plot(0:K-1,ec_nu(1:K,2),0:K-1,scp99_nu(1:K,2),0:K-1,sc2_nu(1:K,2),0:K-1,reg_nu(1:K,2),0:K-1,tracking_nu(1:K,2))
legend('economic','eNMPC-sc w/ \delta=0.99','eNMPC-sc w/ \delta=2.5','regularized','tracking')
xlabel('k')
ylabel('$|x_k - x_{ss}|$')

%avcost=[mean(ec_nu(1:K,3)) , mean(scp99_nu(1:K,3)),mean(sc2_nu(1:K,3)) , mean(reg_nu(1:K,3)) , mean(tracking_nu(1:K,3)) ] 

avprof=-1*[mean(ec_nu(1:K,3)-lecss) , mean(scp99_nu(1:K,3)-lecss),mean(sc2_nu(1:K,3)-lecss) , mean(reg_nu(1:K,3)-lecss) , mean(tracking_nu(1:K,3)-lecss) ] ;

sumcost=[sum(ec_nu(1:K,3)-lecss) , sum(scp99_nu(1:K,3)-lecss),sum(sc2_nu(1:K,3)-lecss) , sum(reg_nu(1:K,3)-lecss) , sum(tracking_nu(1:K,3)-lecss) ] 

figure
bar(avprof)
ylim([ .1 .11])
ylabel('Average profit')
set(gca,'xticklabel',{'econ';'\delta=0.99';'\delta=2.5';'reg';'track'})

ec_nu_unc = importdata('ec_nu_unc.tab');
tracking_nu_unc = importdata('tracking_nu_unc.tab');
reg_nu_unc = importdata('reg_nu_unc.tab');
scp99_nu_unc = importdata('scp99_nu_unc.tab');
sc2_nu_unc = importdata('sc2p5_nu_unc.tab');
figure
plot(0:K-1,ec_nu_unc(1:K,2),0:K-1,scp99_nu_unc(1:K,2),0:K-1,sc2_nu_unc(1:K,2),0:K-1,reg_nu_unc(1:K,2),0:K-1,tracking_nu_unc(1:K,2))
legend('economic','eNMPC-sc w/ \delta=0.99','eNMPC-sc w/ \delta=2.5','regularized','tracking')
xlabel('k')
ylabel('$|x_k - x_{ss}|$')

%avcost_unc=[mean(ec_nu_unc(1:K,3)) , mean(scp99_nu_unc(1:K,3)),mean(sc2_nu_unc(1:K,3))  , mean(reg_nu_unc(1:K,3)) , mean(tracking_nu_unc(1:K,3)) ] 

avprof_unc=-1*[mean(ec_nu_unc(1:K,3)-lecss) , mean(scp99_nu_unc(1:K,3)-lecss),mean(sc2_nu_unc(1:K,3)-lecss)  , mean(reg_nu_unc(1:K,3)-lecss) , mean(tracking_nu_unc(1:K,3)-lecss) ] ;

sumcost_unc=[sum(ec_nu_unc(1:K,3)-lecss) , sum(scp99_nu_unc(1:K,3)-lecss),sum(sc2_nu_unc(1:K,3)-lecss)  , sum(reg_nu_unc(1:K,3)-lecss) , sum(tracking_nu_unc(1:K,3)-lecss) ] 

figure
bar(avprof_unc)
ylim([ .09 .11])
ylabel('Average profit')
set(gca,'xticklabel',{'econ';'\delta=0.99';'\delta=2.5';'reg';'track'})
