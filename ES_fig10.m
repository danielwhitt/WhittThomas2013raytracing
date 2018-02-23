% script to Solve Forced Eliassen-Sawyer Problem in arbitrary background
% flow as in Whitt and Thomas (2013) JPO
% Dan Whitt 

 
 Jplt = 3:3:397;
 Iplt = 3:3:397;

 F2 = F2(Jplt,Iplt);
 S2 = S2(Jplt,Iplt);
 N2 = N2(Jplt,Iplt);
 y_g = y_g(Jplt+1,Iplt+1);
 z_g = z_g(Jplt+1,Iplt+1);
 dz = z_g(2,1)-z_g(1,1)
 dy = y_g(1,2)-y_g(1,1)
 rho0=1025;
 szY = size(F2,2);
 szZ=size(F2,1);
 tau0 = -.24; %N/m^2 % wind stress magnitude
 tau00 = abs(tau0); % this isn't necessarily quite right
 ustar = sqrt(tau00./rho0); % friction velocity
 Hek = .4.*ustar./f; % Turbulent Ekman depth
 DzX = zeros(size(F2));
 DzY = zeros(size(F2));
 nuek = Hek.^2.*f./2; % viscosity in the turbulent Ekman layer
 idxbotEk = ceil(Hek./abs(dz));
 DzX(1:idxbotEk,:) = -tau0./(rho0.*Hek.^2); % vertical derivative of X (momentum source)
 %DzX(10:12,111:113) = -tau0./(rho0.*Hek.^2); % vertical derivative of X (momentum source)

 % add the geostrophic stress 
 %DzX(1:idxbotEk,:) = DzX(1:idxbotEk,:)-nuek.*S2(1:idxbotEk,:)./(f.*Hek.^2);
  %DzY(1:idxbotEk,:) = -tau0./(rho0.*Hek.^2); $ vertical derivative of Y
 %(momentum source)
 
 % change the stratification in the boundary layer 
 %N2OBL = 1e-7;
 %N2(1:(idxbotEk+2),:) = N2OBL;
 %buoy(1:(idxbotEk+2),:) = repmat(buoy((idxbotEk+2),:),[idxbotEk+2 1]) + (z_g(1:(idxbotEk+2),:) - repmat(z_g((idxbotEk+2),:),[(idxbotEk+2) 1])).*N2OBL;
 DyB = zeros(szZ,szY);
 alpha = zeros(szZ,szY);

%nu = .1*(dz^2)*omega/(pi).^2
%nuh = .1*(dy^2)*omega/(pi).^2
nuh = 0;
nuz = 2.0e-5.*ones(szZ,szY); % friction to equilibrate the wavenumber cascade.
%nuz = zeros(size(y_g));
nuzvec = vecES(nuz);
% operators
A = [-1.*ones(szY,1),16.*ones(szY,1),-30.*ones(szY,1),16.*ones(szY,1),-1.*ones(szY,1)];
B = [ones(szY,1),-8.*ones(szY,1),zeros(szY,1),8.*ones(szY,1),-1.*ones(szY,1)];
K2y = (1./(12.*dy.^2)).*spdiags(A,[-2,-1,0,1,2],szY,szY);
K1y = (1./(12.*dy)).*spdiags(B,[-2,-1,0,1,2],szY,szY);
% periodic in y
K2y(1,end) = 16./(12.*dy.^2);
K2y(2,end) = -1./(12.*dy.^2);
K2y(end,1) = 16./(12.*dy.^2);
K2y(end-1,1) = -1./(12.*dy.^2);
K2y(1,end-1) = -1./(12.*dy.^2);
K2y(end,2) = -1./(12.*dy.^2);

K1y(1,end) = -8./(12.*dy);
K1y(2,end) = 1./(12.*dy);
K1y(end,1) = 8./(12.*dy);
K1y(end-1,1) = -1./(12.*dy);
K1y(1,end-1) = 1./(12.*dy);
K1y(end,2) = -1./(12.*dy);

Dy = kron(speye(szZ),K1y);
Dyy = kron(speye(szZ),K2y);
B = fliplr(B);
K1z = (1./(12.*dz)).*spdiags(B,[-2,-1,0,1,2],szZ,szZ);
K2z = (1./(12.*dz.^2)).*spdiags(A,[-2,-1,0,1,2],szZ,szZ);
Dzz = kron(K2z,speye(szY));
Dz = kron(K1z,speye(szY));
N2vec = vecES(N2);
S2vec = vecES(S2);
F2vec = vecES(F2);

Diffz = sparse(Dz*(bsxfun(@times,nuzvec,Dz)));
%RHS of E-S equation
b = vecES(-1i.*omega.*DzY - f.*DzX - DyB - 2.*alpha.*S2) - Diffz*vecES(DzY);
% Elliassen-Sawyer Operator
A = sparse(bsxfun(@times,F2vec,Dzz) + 2.*bsxfun(@times,S2vec,Dz*Dy) + bsxfun(@times,N2vec,Dyy));
% add frequency component (non-hydrostatic)
A  = A - sparse((omega.^2).*Dzz) - sparse((omega.^2).*Dyy);
% add vertical frictional operator to equilibrate the secondary circulation scales
A = A + sparse(2i.*omega.*(Diffz*Dzz)) + sparse(Diffz*Diffz*Dzz) +...
    sparse(2i.*omega.*(Diffz*Dyy)) + sparse(Diffz*Diffz*Dyy) +...
    sparse((1i.*omega).*bsxfun(@times,sparse(Dzz*nuzvec),Dzz)) +...
    sparse((1i.*omega).*bsxfun(@times,sparse(Dz*nuzvec),sparse(Dz*Dzz)))+...
    sparse(Diffz*(sparse(bsxfun(@times,sparse(Dzz*nuzvec),Dzz)))) +...
    sparse(Diffz*(sparse(bsxfun(@times,sparse(Dz*nuzvec),Dz*Dzz))));
%+ sparse(2i.*omega.*nuh.*Dzzyy);
%A = A + 2.*nuh.*nuz.*(Dzzyy*Dzz) + (nuh.*nuh).*(Dyyyy*Dzz);
display('Solve Ax=b exactly')
psivec = A\b;

clear A
psi = matESH(psivec);
%close all

% setup U eqn

dpsidzvec = Dz*psivec;
dpsidyvec = Dy*psivec;
w = dpsidyvec;
v = -dpsidzvec;
%pause
%% plot
%figure
subplot(3,1,3),...
%contourf(y_g./1000+ymean./1000,z_g,real(matESH(w).*86400),100,'linestyle','none'); shading flat
contourf(y_g./1000+ymean./1000,z_g,real(matESH(psivec))./max(abs(psivec(:))),100,'linestyle','none'); 
hold on
caxis(caxis)

hold on
caxis(caxis)
contour(y_g./1000+ymean./1000,z_g,omegamin(Jplt,Iplt),[omega 1e6],'linewidth',2,'color','magenta');
contour(y_g./1000+ymean./1000,z_g,-rhoref.*bg(Jplt+1,Iplt+1)./g -1000,linspace(25.0,25.9,10),'linewidth',1,'color','k');
contour(y_g./1000+ymean./1000,z_g,-rhoref.*bg(Jplt+1,Iplt+1)./g -1000,linspace(26.1,27,10),'linewidth',1,'color','k');
[c1,h1] = contour(y_g./1000+ymean./1000,z_g,-rhoref.*bg(Jplt+1,Iplt+1)./g - 1000,[26.0 1000],'linewidth',2,'color','k');
clabel(c1,h1,'FontSize',12,'FontName','Arial','FontWeight','bold');
ylabel('Depth [m]','FontSize',12,'FontName','Arial','FontWeight','bold');
xlabel('Cross-Stream [km]','FontSize',12,'FontName','Arial','FontWeight','bold');
cbh = colorbar();
set(gca,'FontSize',12,'FontName','Arial','FontWeight','bold');
title('(C) Normalized stream function, \omega = 1.13f','FontSize',12,'FontWeight','bold','FontName','Arial');
xlim([-6.9 16.9])
ylim([-183 -36])
caxis([-.5 .5])