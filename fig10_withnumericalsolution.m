%%  Figure Generation Script for Whitt et al. (2018) 
% Interaction of super-inertial waves with submesoscale cyclonic filaments
% in the North Wall of the Gulf Stream JPO
% JPOD170079
%% WARNING: this is not designed to be a black box. Tuning is required to apply this to other problems.
%% Dan Whitt (dwhitt@ucar.edu) Feb 22 2018
%% Copyright Dan Whitt (dwhitt@ucar.edu) 
%% Written With Matlab v2017b
%% Figure 10
clear all
close all
restoredefaultpath;


g = 9.81; %gravity
rhoref = 1025; % rho_0
fparam = 1.454.*1e-4;
Rearth = 6378.1E3; % meters
%%
load('fig10datafile.mat')

figure;

subplot(3,1,1),...
contourf(y_g(Jplt,Iplt)./1000 +ymean./1000,z_g(Jplt,Iplt),omegamin(Jplt,Iplt)./f,linspace(0,2,61),'linestyle','none');
hold on
caxis([.6 1.6])
colormap(jet)
caxis(caxis)
contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),-rhoref.*bg(Jplt,Iplt)./g -1000,linspace(25.0,25.9,10),'linewidth',1,'color','k');
contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),-rhoref.*bg(Jplt,Iplt)./g -1000,linspace(26.1,27,10),'linewidth',1,'color','k');
contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),omegamin(Jplt,Iplt),[omega 1e6],'linewidth',2,'color','magenta');
contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),omegamin(Jplt,Iplt),[f 1e6],'linewidth',1,'color',[.5 .5 .5]);

[c1,h1] = contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),-rhoref.*bg(Jplt,Iplt)./g - 1000,[26.0 1000],'linewidth',2,'color','k');
clabel(c1,h1,'FontSize',12,'FontName','Arial','FontWeight','bold');
set(gca,'fontsize',12,'fontweight','bold');
ylabel('Depth [m]','FontSize',12,'FontName','Arial','FontWeight','bold');
%xlabel('Cross-Stream [km]','FontSize',20,'FontName','Arial','FontWeight','bold');
cbh = colorbar();
ylabel(cbh,'\omega_{min}','FontSize',12,'FontName','Arial','FontWeight','bold');
title('(A) Minimum Frequency','FontSize',12,'FontWeight','bold','FontName','Arial');

  

subplot(3,1,2),...
contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),-rhoref.*bg(Jplt,Iplt)./g -1000,linspace(25.0,25.9,10),'linewidth',1,'color','k');
hold on
contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),-rhoref.*bg(Jplt,Iplt)./g -1000,linspace(26.1,27,10),'linewidth',1,'color','k');
contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),omegamin(Jplt,Iplt),[omega 1e6],'linewidth',2,'color','magenta');
%contour(y_g(I,J)./1000,z_g(I,J),sqrt(F2)./f,omega_t./f,'linewidth',2,'color','green');
[c1,h1] = contour(y_g(Jplt,Iplt)./1000+ymean./1000,z_g(Jplt,Iplt),-rhoref.*bg(Jplt,Iplt)./g - 1000,[26.0 1000],'linewidth',2,'color','k');
clabel(c1,h1,'FontSize',12,'FontName','Arial','FontWeight','bold');
ylabel('Depth [m]','FontSize',12,'FontName','Arial','FontWeight','bold');
xlabel('Cross-Stream [km]','FontSize',12,'FontName','Arial','FontWeight','bold');
caxis([0 1])
cbh = colorbar();
ylabel(cbh,'|c_g|/|c_g|_{max}','FontSize',12,'FontName','Arial','FontWeight','bold');
set(gca,'FontSize',12,'FontName','Arial','FontWeight','bold');
title('(B) Ray paths, \omega = 1.13f','FontSize',12,'FontWeight','bold','FontName','Arial');

%% Ray tracing
chstart = 1;
dt = 2.5
for i = [-10000,1000]

raystarty = i;
if chstart == 1
    raystartz = -70 -max(z_g(Jplt,10));
else
    raystartz = -70 -max(z_g(Jplt,10));
end
tsteps = 160000;
phasepts = round(1:tstepsperphase:tsteps);

d2udy2 = zeros(398);
d2udz2 = d2udy2;
d2bdy2 = d2udy2;
d2bdzdy = d2udy2;
d2udzdy = d2udy2;
d2bdz2 = d2udy2;

thresh = .1E-2;
try
[y_ray z_ray v_ray w_ray u_ray l_ray m_ray alphaout cgy_ray cgz_ray e_ray s_Mray s_Bray] ...
= raytraceR(z_g(Jplt,Iplt)-max(z_g(Jplt,10)),-dz,y_g(Jplt,Iplt),dy,...
    raystartz,raystarty,Lz,Ly,tsteps,dt,thresh,1/minv,omega, ...
    F2(Jplt,Iplt),S2(Jplt,Iplt),N2(Jplt,Iplt),v0,f,chstart,s_M(Jplt,Iplt),ug(Jplt,Iplt),...
    d2udy2(Jplt,Iplt),d2udz2(Jplt,Iplt),d2bdz2(Jplt,Iplt),d2bdy2(Jplt,Iplt),d2bdzdy(Jplt,Iplt),d2udzdy(Jplt,Iplt));

h10 = scatter(y_ray(1:100:tsteps)./1000+ymean./1000,z_ray(1:100:tsteps)+max(z_g(Jplt,10)),15,...
    sqrt((cgz_ray(1:100:tsteps).^2 + cgy_ray(1:100:tsteps).^2))./max(sqrt((cgz_ray(1:100:tsteps).^2 + cgy_ray(1:100:tsteps).^2))),'filled');
catch
    display('error')
    display(i)
    pause
end

end

pause(1)
chstart = 2;
for i = [-10000,1000]

raystarty = i;
if chstart == 1
    raystartz = -70-max(z_g(Jplt,10));
else
    raystartz = -70-max(z_g(Jplt,10));
end

tsteps = 160000;

d2udy2 = zeros(398);
d2udz2 = d2udy2;
d2bdy2 = d2udy2;
d2bdzdy = d2udy2;
d2udzdy = d2udy2;
d2bdz2 = d2udy2;
thresh = .1E-2;
[y_ray z_ray v_ray w_ray u_ray l_ray m_ray alphaout cgy_ray cgz_ray e_ray s_Mray s_Bray] ...
= raytraceR(z_g(Jplt,Iplt)-max(z_g(Jplt,10)),-dz,y_g(Jplt,Iplt),dy,...
    raystartz,raystarty,Lz,Ly,tsteps,dt,thresh,1/minv,omega, ...
    F2(Jplt,Iplt),S2(Jplt,Iplt),N2(Jplt,Iplt),v0,f,chstart,...
    s_M(Jplt,Iplt),ug(Jplt,Iplt),d2udy2(Jplt,Iplt),d2udz2(Jplt,Iplt),...
    d2bdz2(Jplt,Iplt),d2bdy2(Jplt,Iplt),d2bdzdy(Jplt,Iplt),d2udzdy(Jplt,Iplt));

h10 = scatter(y_ray(1:100:tsteps)./1000+ymean./1000,z_ray(1:100:tsteps)+max(z_g(Jplt,10)),15,...
    sqrt((cgz_ray(1:100:tsteps).^2 + cgy_ray(1:100:tsteps).^2))./max(sqrt((cgz_ray(1:100:tsteps).^2 + cgy_ray(1:100:tsteps).^2))),'filled');
end

%% 


set(gcf,'color','w');
ES_fig10;




