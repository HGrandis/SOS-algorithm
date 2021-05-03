
%  modified SOS algorithm
%  with log resisitivity and linear thickness as parameter
%  and operates in log10 of the resistivity
% (c) Hendra Grandis and Sungkono 2021

clc;

%  CSAMT parameters and other constants
dipl = 2000;        %  Tx dipole length in meters
curr = 10;          %  Tx current in Amperes
xpos = 0;           %  Rx position X in meters
ypos = 5000;        %  Rx position Y in meters
mu = 4*pi*1e-7;
r2d = (180/pi);

%  synthetic data and synthetic model for reference
data = load('data_1.txt');
periods = data(:,1);
rho_obs = data(:,2);
pha_obs = data(:,3);
ndat = length(periods);

model = load('model_1.txt');
res_mod = model(:,1)';
thi_mod = model(:,2)';
lay_mod = length(res_mod);
dmax = 5000;        %  maximum depth for plotting

%  "a priori" resistivity (in log10 Ohm.m) and thickness (in meters) intervals
res_min = 0;
res_max = 3.3;
thi_min = 200;
thi_max = 2000;

%  number of layer, parameter, misfit threshold
nlay = 3;
npar = 2*nlay -1;
rms_max = 0.06;

%  number of models in population and generations or iterations
npop = 100;
ngen = 300;

%  generate initial population randomly
icon = 0;
gre = 10e+10;
for ipop = 1:npop  
    
    for ilay = 1:nlay
        res(ipop,ilay) = res_min + rand*(res_max-res_min);
        thi(ipop,ilay) = thi_min + rand*(thi_max-thi_min);
    end    

%  calculate initial population fitness (misfit)   
    rho_cal = zeros(ndat,1);
    pha_cal = zeros(ndat,1);
    resist = 10.^res(ipop,:);
    thicks = thi(ipop,:);
    
    for idat = 1:ndat
        t = periods(idat);
        zxy = csamt_mex(t,xpos,ypos,dipl,curr,resist,thicks,nlay);
	    zre = real(zxy);
	    zim = imag(zxy);
        w = 2*pi/t;
        rho_cal(idat) = (zre^2+zim^2)/(w*mu);
        pha_cal(idat) = atan2(zim,zre)*r2d; 
    end
    rms(ipop) = misfit(rho_obs,pha_obs,rho_cal,pha_cal);
    
%  find best model in initial population
    if (rms(ipop) < gre)
       gre = rms(ipop);
       res_best = resist;
       thi_best = thicks;
       rms_best = rms(ipop);
    end
    
%  isolate model with misfit lower than threshold value    
    if (rms(ipop) < rms_max)
        icon = icon + 1;
        res_fit(icon,:) = res(ipop,:);
        thi_fit(icon,:) = thi(ipop,:);
    end
end
%  end of initial population generation

%  start SOS algorithm
%  loop over generations
for igen = 1:ngen
   
    for i = npop
        
%  mutualism phase
%  select two models other than i-th model
        do = true; 
        while do
           j = randi(npop);
           if  (j ~= i), break, end 
        end
        do = true; 
        while do
           k = randi(npop);
           if (k ~= i) && (k ~= j), break, end 
        end
%  save misfit of (old) models i and j
        rms_old1 = rms(i);
        rms_old2 = rms(j);
                
%  mutual vector for resistivity and thickness of models i and j
        mv_res = 0.5*(res(i,:) + res(j,:)); 
        mv_thi = 0.5*(thi(i,:) + thi(j,:));
          
%  perturbation of i-th model relative to k-th model
        res_new1 = res(i,:) + rand(1,nlay).*(res(k,:) - mv_res);
        thi_new1 = thi(i,:) + rand(1,nlay).*(thi(k,:) - mv_thi);
        
%  limit perturbation to "a priori" parameters' intervals
        for ilay = 1:nlay
            if (res_new1(ilay) < res_min) 
                res_new1(ilay) = res_min;
            end
            if (res_new1(ilay) > res_max) 
                res_new1(ilay) = res_max;
            end  
            if (thi_new1(ilay) < thi_min) 
                thi_new1(ilay) = thi_min;
            end
            if (thi_new1(ilay) > thi_max) 
                thi_new1(ilay) = thi_max;
            end
        end
        
%  calculate misfit of new i-th model 
        resist = 10.^res_new1;
        thicks = thi_new1;
        for idat = 1:ndat
            t = periods(idat);
            zxy = csamt_mex(t,xpos,ypos,dipl,curr,resist,thicks,nlay);
	        zre = real(zxy);
	        zim = imag(zxy);
            w = 2*pi/t;
            rho_cal(idat) = (zre^2+zim^2)/(w*mu);
            pha_cal(idat) = atan2(zim,zre)*r2d;             
        end
        rms_new1 = misfit(rho_obs,pha_obs,rho_cal,pha_cal);
        
%  replace i-th old model with new one if it has lower misfit
        if (rms_new1 < rms_old1)
            res(i,:) = res_new1;
            thi(i,:) = thi_new1;
            rms(i) = rms_new1;
        end
           
%  perturbation of j-th model relative to k-th model    
        res_new2 = res(j,:) + rand(1,nlay).*(res(k,:) - mv_res);
        thi_new2 = thi(j,:) + rand(1,nlay).*(thi(k,:) - mv_thi);

%  limit perturbation to "a priori" parameters' intervals        
        for ilay = 1:nlay
            if (res_new2(ilay) < res_min) 
                res_new2(ilay) = res_min;
            end
            if (res_new2(ilay) > res_max) 
                res_new2(ilay) = res_max;
            end            
            if (thi_new2(ilay) < thi_min)
                thi_new2(ilay) = thi_min;
            end
            if (thi_new2(ilay) > thi_max) 
                thi_new2(ilay) = thi_max;
            end
        end
        
%  calculate misfit of new j-th model
        resist = 10.^res_new2;
        thicks = thi_new2;
        for idat = 1:ndat
            t = periods(idat);
            zxy = csamt_mex(t,xpos,ypos,dipl,curr,resist,thicks,nlay);
	        zre = real(zxy);
	        zim = imag(zxy);
            w = 2*pi/t;
            rho_cal(idat) = (zre^2+zim^2)/(w*mu);
            pha_cal(idat) = atan2(zim,zre)*r2d;
        end
        rms_new2 = misfit(rho_obs,pha_obs,rho_cal,pha_cal);
        
%  replace j-th old model with new one if it has lower misfit
        if (rms_new2 < rms_old2)
            res(j,:) = res_new2;
            thi(j,:) = thi_new2;
            rms(j) = rms_new2;
        end     
%  end of mutualism phase

%  compare the new models with the best model 
        if (rms(i) < rms_best)
           res_best = res(i,:);
           thi_best = thi(i,:);
           rms_best = rms(i);
        end
        if (rms(j) < rms_best)
           res_best = res(j,:);
           thi_best = thi(j,:);
           rms_best = rms(j);
        end
        
%  commensalism phase
%  select j-th model other than i-th model
        do = true; 
        while do
           j = randi(npop);
           if  (j ~= i), break, end 
        end
%  save misfit of (old) models i
        rms_old1 = rms(i);
        
%  perturbation of i-th model relative to best and j-th model
        rand1 = -0.5 + rand(1,nlay);
        rand2 = -0.5 + rand(1,nlay);
%         rand1 = 0.4 + 0.5*rand(1,nlay);
%         rand2 = 0.4 + 0.5*rand(1,nlay);
        res_new1 = res(i,:) + rand1.*(res_best - res(j,:));
        thi_new1 = thi(i,:) + rand2.*(thi_best - thi(j,:));

%  limit perturbation to "a priori" parameters' intervals
        for ilay = 1:nlay
            if (res_new1(ilay) < res_min) 
                res_new1(ilay) = res_min;
            end
            if (res_new1(ilay) > res_max) 
                res_new1(ilay) = res_max;
            end            
            if (thi_new1(ilay) < thi_min)
                thi_new1(ilay) = thi_min;
            end
            if (thi_new1(ilay) > thi_max) 
                thi_new1(ilay) = thi_max;
            end
        end   
        
%  calculate misfit of new i-th model 
        resist = 10.^res_new1;
        thicks = thi_new1;
        for idat = 1:ndat
            t = periods(idat);
            zxy = csamt_mex(t,xpos,ypos,dipl,curr,resist,thicks,nlay);
	        zre = real(zxy);
	        zim = imag(zxy);
            w = 2*pi/t;
            rho_cal(idat) = (zre^2+zim^2)/(w*mu);
            pha_cal(idat) = atan2(zim,zre)*r2d;
        end
        rms_new1 = misfit(rho_obs,pha_obs,rho_cal,pha_cal);
        
%  replace old i-th model with new one if it has lower misfit
        if (rms_new1 < rms_old1)
            res(i,:) = res_new1;
            thi(i,:) = thi_new1;
            rms(i) = rms_new1;
        end
%  end commensalism phase
        
%  parasitism phase
%  select model other than i-th model
        do = true; 
        while do
           j = randi(npop);
           if  (j ~= i), break, end 
        end
%  save misfit of (old) models i and j
        rms_old1 = rms(i);
        rms_old2 = rms(j);
        
%  perturbation of i-th model as parasite
        res_new1 = res(i,:);
        thi_new1 = thi(i,:);
        rand1 = randi(npar);
        if (rand1 <= nlay)
            res_new1(rand1) = res_min + rand*(res_max-res_min);
        else
            ipar = rand1 - nlay;
            thi_new1(ipar) = thi_min + rand*(thi_max-thi_min);
        end
        
%  calculate misfit of new i-th model as parasite
        resist = 10.^res_new1;
        thicks = thi_new1;
        for idat = 1:ndat
            t = periods(idat);
            zxy = csamt_mex(t,xpos,ypos,dipl,curr,resist,thicks,nlay);
	        zre = real(zxy);
	        zim = imag(zxy);
            w = 2*pi/t;
            rho_cal(idat) = (zre^2+zim^2)/(w*mu);
            pha_cal(idat) = atan2(zim,zre)*r2d;
        end
        rms_new1 = misfit(rho_obs,pha_obs,rho_cal,pha_cal);
        
%  replace j-th model with parasite if it has lower misfit
        if (rms_new1 < rms_old2)
            res(j,:) = res_new1;
            thi(j,:) = thi_new1;
            rms(j) = rms_new1;
        end
%  end of parasitism phase

    end
%  end of loop over population in a generation / iteration

%  find best model in population at igen-th iteration
    gre = 10e+10;
    for i = 1:npop 
        if (rms(i) < gre)
           gre = rms(i);
           res_best = res(i,:);
           thi_best = thi(i,:);
           rms_best = rms(i);
        end
        
%  isolate model with misfit lower than threshold value 
        if (rms(i) < rms_max)
           icon = icon + 1;
           res_fit(icon,:) = res(i,:);
           thi_fit(icon,:) = thi(i,:);
        end
    end
    
%  calculate best model response
    resist = 10.^res_best;
    thicks = thi_best;
    for idat = 1:ndat
        t = periods(idat);
        zxy = csamt_mex(t,xpos,ypos,dipl,curr,resist,thicks,nlay);
	    zre = real(zxy);
	    zim = imag(zxy);
        w = 2*pi/t;
        rho_cal(idat) = (zre^2+zim^2)/(w*mu);
        pha_cal(idat) = atan2(zim,zre)*r2d;
    end
  
    rms_iter(igen) = misfit(rho_obs,pha_obs,rho_cal,pha_cal);
    iter(igen) = igen;
    
%  plot best model at each iteration    
    figure(1)
    set(figure(1),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gcf,'color','w');

    subplot(2,2,1);
    loglog(periods,rho_obs,'.r','markers',15);
    hold on
    loglog(periods,rho_cal,'b','LineWidth',2);

    set(gca,'FontSize',11,'LineWidth',1)
    xlabel('Period (sec.)','fontweight','bold','fontsize',12);
    ylabel('App. Resistivity (Ohm.m)','fontweight','bold','fontsize',12) ;
    axis([0.001 10 1 10000])
    hold off

    subplot(2,2,3);
    semilogx(periods,pha_obs,'.r','markers',15);
    hold on
    semilogx(periods,pha_cal,'b','LineWidth', 2);

    set(gca,'FontSize',11,'LineWidth',1)
    xlabel('Period (sec.)','fontweight','bold','fontsize',12);
    ylabel('Phase (deg.)','fontweight','bold','fontsize',12) ;
    axis([0.001 10 0 90])
    hold off

    rr = [0,10.^res_best];
    tt = [0,cumsum(thi_best(1:nlay-1)),dmax]; 
    rrr = [0,res_mod];
    ttt = [0,cumsum(thi_mod(1:lay_mod-1)),dmax]; 
    subplot(2,2,[2,4]), 
    stairs(rr,tt,'b-','LineWidth', 2);
    hold on
    stairs(rrr,ttt,'r--','LineWidth', 1);
    set(gca,'Ydir','reverse'); 
    set(gca,'Xscale','log'); 
    axis([0.1 10000 0 dmax]); 
    set(gca,'FontSize',11,'LineWidth',1)
    xlabel('Resistivity (Ohm.m)','fontweight','bold','fontsize',12); 
    ylabel('Depth (m)','fontweight','bold','fontsize',12);
    hold off
    drawnow
end
%  end of iterations over generations

%  calculate average model with misfit lower than threshold value
if icon ~= 0
   for ilay = 1:nlay
       res_aver(ilay) = mean(res_fit(:,ilay));
       thi_aver(ilay) = mean(thi_fit(:,ilay));
   end
else
%  use best model otherwise    
   res_aver = res_best;
   thi_aver = thi_best;
end

%  calculate misfit of average or best model 
resist = 10.^res_aver;
thicks = thi_aver;
for idat = 1:ndat
    t = periods(idat);
    zxy = csamt_mex(t,xpos,ypos,dipl,curr,resist,thicks,nlay);
	zre = real(zxy);
	zim = imag(zxy);
    w = 2*pi/t;
    rho_cal(idat) = (zre^2+zim^2)/(w*mu);
    pha_cal(idat) = atan2(zim,zre)*r2d;
end
rms_aver = misfit(rho_obs,pha_obs,rho_cal,pha_cal)
num_aver = icon

figure(2)
set(figure(2),'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'color','w');

subplot(2,2,1);
loglog(periods,rho_obs,'.r','markers',15);
hold on
loglog(periods,rho_cal,'b','LineWidth', 2);

set(gca,'FontSize',11,'LineWidth',1)

xlabel('Period (sec.)','fontweight','bold','fontsize',12);
ylabel('App. Resistivity (Ohm.m)','fontweight','bold','fontsize',12) ;
axis([0.001 10 1 10000])

subplot(2,2,3);
semilogx(periods,pha_obs,'.r','markers',15);
hold on
semilogx(periods,pha_cal,'b','LineWidth', 2);

set(gca,'FontSize',11,'LineWidth',1)

xlabel('Period (sec.)','fontweight','bold','fontsize',12);
ylabel('Phase (deg.)','fontweight','bold','fontsize',12) ;
axis([0.001 10 0 90])

rr = [0,10.^res_aver];
tt = [0,cumsum(thi_aver(1:nlay-1)),dmax];
rrr = [0,res_mod];
ttt = [0,cumsum(thi_mod(1:lay_mod-1)),dmax]; 
subplot(2,2,[2,4]), 
stairs(rr,tt,'b-','LineWidth', 2);
hold on
stairs(rrr,ttt,'r--','LineWidth', 1);

set(gca,'Ydir','reverse'); 
set(gca,'Xscale','log'); 
axis([0.1 10000 0 dmax]); 

set(gca,'FontSize',11,'LineWidth',1)
xlabel('Resistivity (Ohm.m)','fontweight','bold','fontsize',12); 
ylabel('Depth (m)','fontweight','bold','fontsize',12);
hold on

figure(3)
plot(iter,rms_iter,'b','LineWidth', 2);
xlabel('iterations','fontweight','bold','fontsize',12);
ylabel('misfit','fontweight','bold','fontsize',12) ;
axis([1 ngen 0 0.5])

%  Save results to txt file
%  Calculated / theoretical response of inverse model
file_dat = fopen('datcal_1.txt','w');
for i = 1:ndat
    fprintf(file_dat,'%12.5f %12.4f %12.4f \n',periods(i),rho_cal(i),pha_cal(i));
end
fclose(file_dat);

%  RMS error as function of iteratiions for plotting
file_dat = fopen('rmsite_1.txt','w');
for i = 1:ngen
    fprintf(file_dat,'%12.1f %12.4f \n',iter(i),rms_iter(i));
end
fclose(file_dat);

%  Final RMS error and number of averaged models
file_dat = fopen('rmsfin_1.txt','w');
fprintf(file_dat,'%12.1f %12.4f \n',num_aver,rms_aver);
fclose(file_dat);

%  Inverse model for plotting
file_dat = fopen('modinv_1.txt','w');
fprintf(file_dat,'%12.1f %12.1f \n',resist(1),0);
depth = 0;
for i = 1:nlay-1
    depth = depth + thi_aver(i);
    fprintf(file_dat,'%12.1f %12.4f \n',10.^res_aver(i),depth/1000);
    fprintf(file_dat,'%12.1f %12.4f \n',10.^res_aver(i+1),depth/1000);
end
fprintf(file_dat,'%12.1f %12.1f \n',10.^res_aver(nlay),dmax/1000);
fclose(file_dat);

%  All models below threshold RMS error for plotting
file_dat = fopen('modall_1.txt','w');
for j = 1:icon
fprintf(file_dat,'%12.1f %12.1f \n',10^res_fit(j,1),0);
depth = 0;
for i = 1:nlay-1
    depth = depth + thi_fit(j,i);
    fprintf(file_dat,'%12.1f %12.4f \n',10^res_fit(j,i),depth/1000);
    fprintf(file_dat,'%12.1f %12.4f \n',10^res_fit(j,i+1),depth/1000);
end
fprintf(file_dat,'%12.1f %12.4f \n',10^res_fit(j,nlay),dmax/1000);
fprintf(file_dat,'%s \n', '>>');
end
fclose(file_dat);