function [y_ray z_ray v_ray w_ray u_ray l_ray m_ray alphaout cgy_ray cgz_ray e_ray s_Mray s_Bray] = raytraceR(Z,dz,Y,dy,...
   startz,starty,Lz,Ly,tsteps,dt,thresh,m0,omega, ...
    F2,S2,N2,v0,f,chstart,s_M,ug,d2udy2,d2udz2,d2bdz2,d2bdy2,d2bdzdy,d2udzdy)
%%% GENERATION PAPER VERSION %%%
%%%% REFLECTION RULES:%%%%%%% 
% Internal reflection: switch characteristic and hold the sign of m fixed
% Bottom reflection: stop ray, do not reflect
% Top reflection: switch characteristic, hold sign of m fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function creates y and z coordinates of rays.
% Input: Z coordinate (matrix)
% Input: Y coordinate (matrix)
% 
% Input: W group velocity in Z direction
% Input: V group velocity in Y direction
% Input: startz (coordinate location where the ray should start in z) 
% Input: starty (coordinate location where the ray should start in y)
% assumes a grid [-1 0] in Z and a grid of [-1 1] in Y with Zstar and Ystar
% corresponding
% Dan Whitt
% January 10, 2011
% For use by characteristics.m script
% Update As of January 11, 2011 we are computing
% ray velocities v_ag and w with this script as well
% from the polarization relations and returning them as v_ray
% and w_ray

%addpath('/Users/dwhitt/Google Drive/copymatlabstuff/matlabstuff/')

% initialize flags
reflectflag = 0;
%reflectI = 0;
%stopB = 0;
test1flag = 1; % set to 1 to hold m fixed
%anomflag = omega - sqrt(F2); % anomalously low frequency check
om = sqrt(F2 - (S2.^2)./N2);
%trapflag = omega - om;

% initialize ray parameter vectors
y_ray = zeros(tsteps,1);
z_ray = zeros(tsteps,1);
v_ray = zeros(tsteps,1);
w_ray = zeros(tsteps,1);
u_ray = zeros(tsteps,1);
l_ray = zeros(tsteps,1);
l1_ray = zeros(tsteps,1);
l2_ray = zeros(tsteps,1);
m_ray = zeros(tsteps,1);
m1_ray = zeros(tsteps,1);
m2_ray = zeros(tsteps,1);

alphaout = zeros(tsteps,1);
e_ray = zeros(tsteps,1);
cgy_ray = zeros(tsteps,1);
cgz_ray = zeros(tsteps,1);
x1_ray = zeros(tsteps,1); % aspect ratios
x2_ray = zeros(tsteps,1);
ch_ray = zeros(tsteps,1);
cellz = zeros(tsteps,1);
celly = zeros(tsteps,1);
s_Mray = zeros(tsteps,1);
s_Bray = zeros(tsteps,1);
domegady1_ray = zeros(tsteps,1);
domegadz1_ray = zeros(tsteps,1);
domegady2_ray = zeros(tsteps,1);
domegadz2_ray = zeros(tsteps,1);

% start ray trace

% initialize coordinates and start location
ch_ray(1) = chstart;
z_ray(1) = startz;
y_ray(1) = starty;
z2 = zeros(size(Z,1),1);
y2 = zeros(size(Y,2),1);
z2(:) = Z(:,1);
y2(:) = Y(1,:);
szZ = size(Z,1);
szY = size(Y,2);

% interpolate background parameters to ray starting location
F2_bilin = bilininterp(z2,y2,F2,startz,starty);
S2_bilin = bilininterp(z2,y2,S2,startz,starty);
N2_bilin = bilininterp(z2,y2,N2,startz,starty);
% 
d2udy2_bilin = bilininterp(z2,y2,d2udy2,startz,starty);
d2udz2_bilin = bilininterp(z2,y2,d2udz2,startz,starty);
d2bdy2_bilin = bilininterp(z2,y2,d2bdy2,startz,starty);
d2bdz2_bilin = bilininterp(z2,y2,d2bdz2,startz,starty);
d2udzdy_bilin = bilininterp(z2,y2,d2udzdy,startz,starty);
d2bdzdy_bilin = bilininterp(z2,y2,d2bdzdy,startz,starty);
%
s_Bray(1) = S2_bilin./N2_bilin;
s_Mray(1) = bilininterp(z2,y2,s_M,startz,starty);

% compute possible aspect ratios - (nonhydrostatic)
x1_ray(1) = (-S2_bilin + sqrt(S2_bilin.^2 + (omega.^2 - F2_bilin).*(N2_bilin - omega.^2)))./(N2_bilin - omega.^2);
x2_ray(1) = (-S2_bilin - sqrt(S2_bilin.^2 + (omega.^2 - F2_bilin).*(N2_bilin - omega.^2)))./(N2_bilin - omega.^2);

% compute changes in wave number
domegady1_ray(1) = (-f./(2.*omega)).*(d2udy2_bilin - 2.*d2udzdy_bilin.*x1_ray(1) + d2udz2_bilin.*(x1_ray(1)).^2);
domegadz1_ray(1) =  (1./(2.*omega)).*(d2bdy2_bilin  - 2.*d2bdzdy_bilin.*x1_ray(1) + d2bdz2_bilin.*(x1_ray(1)).^2);
domegady2_ray(1) = (-f./(2.*omega)).*(d2udy2_bilin - 2.*d2udzdy_bilin.*x2_ray(1) + d2udz2_bilin.*(x2_ray(1)).^2);
domegadz2_ray(1) =  (1./(2.*omega)).*(d2bdy2_bilin  - 2.*d2bdzdy_bilin.*x2_ray(1) + d2bdz2_bilin.*(x2_ray(1)).^2);

% initial wave numbers
m_ray(1) = m0;
m1_ray(1) = m0;
m2_ray(1) = m0;
l1_ray(1) = x1_ray(1).*m0;
l2_ray(1) = x2_ray(1).*m0;

% hydrostatic group velocities
cgz_1 = (-S2_bilin*l1_ray(1)/(m1_ray(1)^2) - N2_bilin*(l1_ray(1)^2)/(m1_ray(1)^3))/(omega);
cgz_2 = (-S2_bilin*l2_ray(1)/(m2_ray(1)^2) - N2_bilin*(l2_ray(1)^2)/(m2_ray(1)^3))/(omega);
cgy_1 = (S2_bilin/m1_ray(1) + N2_bilin*l1_ray(1)/(m1_ray(1)^2))/(omega);
cgy_2 = (S2_bilin/m2_ray(1) + N2_bilin*l2_ray(1)/(m2_ray(1)^2))/(omega);

% non-hydrostatic group velocities - outdated
%cgz_1n = ((F2_bilin*m_ray(1)*l1_ray(1)^2) + (S2_bilin*l1_ray(1)^3) - ...
%    (S2_bilin*l1_ray(1)*m_ray(1)^2) - N2_bilin*m_ray(1)*l1_ray(1)^2)/(omega*((m_ray(1)^2 + l1_ray(1)^2)^2));
%cgz_2n = ((F2_bilin*m_ray(1)*l2_ray(1)^2) + (S2_bilin*l2_ray(1)^3) - ...
%    (S2_bilin*l2_ray(1)*m_ray(1)^2) - N2_bilin*m_ray(1)*l2_ray(1)^2)/(omega*((m_ray(1)^2 + l2_ray(1)^2)^2));
%cgy_1n = (-F2_bilin*l1_ray(1)*m_ray(1)^2 + S2_bilin*m_ray(1)^3 - S2_bilin*m_ray(1)*l1_ray(1)^2 + ...
 %   N2_bilin*l1_ray(1)*m_ray(1)^2)/(omega*((m_ray(1)^2 + l1_ray(1)^2)^2));
%cgy_2n = (-F2_bilin*l2_ray(1)*m_ray(1)^2 + S2_bilin*m_ray(1)^3 - S2_bilin*m_ray(1)*l2_ray(1)^2 + ...
%    N2_bilin*l2_ray(1)*m_ray(1)^2)/(omega*((m_ray(1)^2 + l2_ray(1)^2)^2));

if ch_ray(1) == 1
    cgz_ray(1) = cgz_1;
    l_ray(1) = l1_ray(1);
    cgy_ray(1) = cgy_1;
    alphaout(1) = x1_ray(1);
elseif ch_ray(1) == 2
    cgz_ray(1) = cgz_2;
    l_ray(1) = l2_ray(1);
    cgy_ray(1) = cgy_2;
    alphaout(1) = x2_ray(1);
end

% Energy Considerations - assuming equipartition
% e_ray(1) = 1;
% E = @(vs) 0.5.*1025.*(vs.^2 + (-l_ray(1).*vs./m_ray(1)).^2) - e_ray(1);
% v0 = fzero(E,v0);


% polarization relations - see written notes
% v_ray(1) = v0;
% w_ray(1) = -l_ray(1).*v_ray(1)./m_ray(1);    % from continuity
% u_ray(1) = v_ray(1).*(F2_bilin + S2_bilin.*l_ray(1)./m_ray(1))./(-1i.*f.*omega); 
%

if imag(sqrt(cgy_ray(1).^2 + cgz_ray(1).^2)) < 1E-9

    delay = 0;
    
    for t = 2:tsteps

        cellz(t-1) = (z_ray(t-1) + Lz)/dz;
        celly(t-1) = (y_ray(t-1) + Ly./2)/dy;
        display('t, zcell, ycell:')
        display(t)
        display(cellz(t-1))
        display(celly(t-1))
        
        if celly(t-1) <= 5 || celly(t-1) >= (szY-5) || cellz(t-1) <= 5
            % we hit the bottom or the sides, BREAK
            masknan = y_ray == 0;
            y_ray(masknan) = nan;
            z_ray(masknan) = nan;
            break
        end
        
        
        if cellz(t-1) > 5 && celly(t-1) > 5 && (cellz(t-1) < (szZ-5)) && celly(t-1) < (szY-5) 
            % Ray is not near a domain  boundary
            
            if (delay > 75 && sqrt(real(cgy_ray(t-1)).^2 + real(cgz_ray(t-1))^2) < 2*thresh) 
                % Should we internally reflect??
                
                m_ray(t) = m_ray(t-5);

                if ch_ray(t-2) == 1
                    ch_ray(t) = 2;
                else
                    ch_ray(t) = 1;
                end
                
                reflectflag = 1
                delay = 0;
                y_ray(t) = y_ray(t-5);
                z_ray(t) = z_ray(t-5);
                
            else
                ch_ray(t) = ch_ray(t-1);
                delay = delay + 1;
                y_ray(t) = y_ray(t-1) + real(cgy_ray(t-1))*dt;
                z_ray(t) = z_ray(t-1) + real(cgz_ray(t-1))*dt;
            end
            
            F2_bilin = bilininterp(z2,y2,F2,z_ray(t),y_ray(t));
            S2_bilin = bilininterp(z2,y2,S2,z_ray(t),y_ray(t));
            N2_bilin = bilininterp(z2,y2,N2,z_ray(t),y_ray(t));
            d2udy2_bilin = bilininterp(z2,y2,d2udy2,z_ray(t),y_ray(t));
            d2udz2_bilin = bilininterp(z2,y2,d2udz2,z_ray(t),y_ray(t));
            d2bdy2_bilin = bilininterp(z2,y2,d2bdy2,z_ray(t),y_ray(t));
            d2bdz2_bilin = bilininterp(z2,y2,d2bdz2,z_ray(t),y_ray(t));
            d2udzdy_bilin = bilininterp(z2,y2,d2udzdy,z_ray(t),y_ray(t));
            d2bdzdy_bilin = bilininterp(z2,y2,d2bdzdy,z_ray(t),y_ray(t));

            s_Mray(t) = bilininterp(z2,y2,s_M,z_ray(t),y_ray(t));
            s_Bray(t) = S2_bilin./N2_bilin;
            
            x1_ray(t) = (-S2_bilin + sqrt(S2_bilin.^2 + (omega.^2 - F2_bilin).*(N2_bilin - omega.^2)))./(N2_bilin - omega.^2);
            x2_ray(t) = (-S2_bilin - sqrt(S2_bilin.^2 + (omega.^2 - F2_bilin).*(N2_bilin - omega.^2)))./(N2_bilin - omega.^2);
            
            domegady1_ray(t) = (-f./(2.*omega)).*(d2udy2_bilin - 2.*d2udzdy_bilin.*x1_ray(t) + d2udz2_bilin.*(x1_ray(t)).^2);
            domegadz1_ray(t) =  (1./(2.*omega)).*(d2bdy2_bilin  - 2.*d2bdzdy_bilin.*x1_ray(t) + d2bdz2_bilin.*(x1_ray(t)).^2);
            domegady2_ray(t) = (-f./(2.*omega)).*(d2udy2_bilin - 2.*d2udzdy_bilin.*x2_ray(t) + d2udz2_bilin.*(x2_ray(t)).^2);
            domegadz2_ray(t) =  (1./(2.*omega)).*(d2bdy2_bilin  - 2.*d2bdzdy_bilin.*x2_ray(t) + d2bdz2_bilin.*(x2_ray(t)).^2);
            
            if test1flag == 1
                if reflectflag == 1
                    m1_ray(t) = m_ray(t);
                    m2_ray(t) = m_ray(t);
                    reflectflag = 0;
                else
                    m1_ray(t) = m_ray(t-1);
                    m2_ray(t) = m_ray(t-1);
                end
                l1_ray(t) = m1_ray(t).*real(x1_ray(t));
                l2_ray(t) = m2_ray(t).*real(x2_ray(t));
            else
                if reflectflag == 0
                    m1_ray(t) = m_ray(t-1) - dt.*real(domegadz1_ray(t));
                    m2_ray(t) = m_ray(t-1) - dt.*real(domegadz2_ray(t));
                    %if real(x1_ray(t)) == 0 
                     %   m1_ray(t) = m1_ray(t-1);
                     %   m2_ray(t) = m2_ray(t-1);
                    %else
                    l1_ray(t) = m1_ray(t).*real(x1_ray(t));
                    l2_ray(t) = m2_ray(t).*real(x2_ray(t));
                    %end
                else
                    m1_ray(t) = m_ray(t-5);
                    m2_ray(t) = m_ray(t-5);
                    l1_ray(t) = m1_ray(t).*real(x1_ray(t-5));
                    l2_ray(t) = m2_ray(t).*real(x2_ray(t-5));
                    reflectflag = 0;
                end
            end
            cgz_1 = (-S2_bilin*l1_ray(t)./(m1_ray(t).^2) - N2_bilin*(l1_ray(t)^2)./(m1_ray(t).^3))/(omega);
            cgz_2 = (-S2_bilin*l2_ray(t)./(m2_ray(t).^2) - N2_bilin*(l2_ray(t)^2)./(m2_ray(t).^3))/(omega);
            cgy_1 = (S2_bilin./m1_ray(t) + N2_bilin*l1_ray(t)./(m1_ray(t)^2))/(omega);
            cgy_2 = (S2_bilin./m2_ray(t) + N2_bilin*l2_ray(t)./(m2_ray(t)^2))/(omega);

% Non-hydrostatic group velocities            
%             cgz_1 = ((F2_bilin*m_ray(t)*l1_ray(t)^2) + (S2_bilin*l1_ray(t)^3) - ...
%                 (S2_bilin*l1_ray(t)*m_ray(t)^2) - N2_bilin*m_ray(t)*l1_ray(t)^2)/(omega*((m_ray(t)^2 + l1_ray(t)^2)^2));
%             cgz_2 = ((F2_bilin*m_ray(t)*l2_ray(t)^2) + (S2_bilin*l2_ray(t)^3) - ...
%                 (S2_bilin*l2_ray(t)*m_ray(t)^2) - N2_bilin*m_ray(t)*l2_ray(t)^2)/(omega*((m_ray(t)^2 + l2_ray(t)^2)^2));
%             cgy_1 = (-F2_bilin*l1_ray(t)*m_ray(t)^2 + S2_bilin*m_ray(t)^3 - S2_bilin*m_ray(t)*l1_ray(t)^2 + ...
%                 N2_bilin*l1_ray(t)*m_ray(t)^2)/(omega*((m_ray(t)^2 + l1_ray(t)^2)^2));
%             cgy_2 = (-F2_bilin*l2_ray(t)*m_ray(t)^2 + S2_bilin*m_ray(t)^3 - S2_bilin*m_ray(t)*l2_ray(t)^2 + ...
%                 N2_bilin*l2_ray(t)*m_ray(t)^2)/(omega*((m_ray(t)^2 + l2_ray(t)^2)^2));
            
            if ch_ray(t) == 1
                cgz_ray(t) = cgz_1;
                l_ray(t) = l1_ray(t);
                m_ray(t) = m1_ray(t);
                cgy_ray(t) = cgy_1;
                alphaout(t) = x1_ray(t);
            elseif ch_ray(t) == 2
                cgz_ray(t) = cgz_2;
                l_ray(t) = l2_ray(t);
                m_ray(t) = m2_ray(t);
                cgy_ray(t) = cgy_2;
                alphaout(t) = x2_ray(t);
            end
            
            % Energy Considerations
%             e_ray(t) = e_ray(1).*sqrt(cgy_ray(1).^2 + cgz_ray(1).^2)./sqrt(cgz_ray(t).^2 + cgy_ray(t).^2);
%             E = @(vs) 0.5.*1025.*(vs.^2 + (-l_ray(t).*vs./m_ray(t)).^2) - e_ray(t);
%             v0 = fzero(E,v_ray(t-1));
            
            % polarization relations
            v_ray(t) = v0;
            w_ray(t) = -l_ray(t).*v_ray(t)./m_ray(t);    % from continuity
            u_ray(t) = v_ray(t).*(F2_bilin + S2_bilin.*l_ray(t)./m_ray(t))./(-1i.*f.*omega);  
            
        else
            %if flag(round(cellz(t-3)),round(celly(t-3))) > 0

            
            m_ray(t) = m_ray(t-5);
            
            
            if ch_ray(t-2) == 1
                ch_ray(t) = 2;
            else
                ch_ray(t) = 1;
            end
            reflectflag = 1
            delay = 0;

            y_ray(t) = y_ray(t-5);
            z_ray(t) = z_ray(t-5);
            F2_bilin = bilininterp(z2,y2,F2,z_ray(t),y_ray(t));
            S2_bilin = bilininterp(z2,y2,S2,z_ray(t),y_ray(t));
            N2_bilin = bilininterp(z2,y2,N2,z_ray(t),y_ray(t));
            
            d2udy2_bilin = bilininterp(z2,y2,d2udy2,z_ray(t),y_ray(t));
            d2udz2_bilin = bilininterp(z2,y2,d2udz2,z_ray(t),y_ray(t));
            d2bdy2_bilin = bilininterp(z2,y2,d2bdy2,z_ray(t),y_ray(t));
            d2bdz2_bilin = bilininterp(z2,y2,d2bdz2,z_ray(t),y_ray(t));
            d2udzdy_bilin = bilininterp(z2,y2,d2udzdy,z_ray(t),y_ray(t));
            d2bdzdy_bilin = bilininterp(z2,y2,d2bdzdy,z_ray(t),y_ray(t));
            
            s_Mray(t) = bilininterp(z2,y2,s_M,z_ray(t),y_ray(t));
            s_Bray(t) = S2_bilin./N2_bilin;
            
            x1_ray(t) = (-S2_bilin + sqrt(S2_bilin.^2 + (omega.^2 - F2_bilin).*(N2_bilin - omega.^2)))./(N2_bilin - omega.^2);
            x2_ray(t) = (-S2_bilin - sqrt(S2_bilin.^2 + (omega.^2 - F2_bilin).*(N2_bilin - omega.^2)))./(N2_bilin - omega.^2);
            
            domegady1_ray(t) = (-f./(2.*omega)).*(d2udy2_bilin - 2.*d2udzdy_bilin.*x1_ray(t) + d2udz2_bilin.*(x1_ray(t)).^2);
            domegadz1_ray(t) =  (1./(2.*omega)).*(d2bdy2_bilin  - 2.*d2bdzdy_bilin.*x1_ray(t) + d2bdz2_bilin.*(x1_ray(t)).^2);
            domegady2_ray(t) = (-f./(2.*omega)).*(d2udy2_bilin - 2.*d2udzdy_bilin.*x2_ray(t) + d2udz2_bilin.*(x2_ray(t)).^2);
            domegadz2_ray(t) =  (1./(2.*omega)).*(d2bdy2_bilin  - 2.*d2bdzdy_bilin.*x2_ray(t) + d2bdz2_bilin.*(x2_ray(t)).^2);
            
           
            if test1flag == 1
                if reflectflag == 1
                    m1_ray(t) = m_ray(t);
                    m2_ray(t) = m_ray(t);
                    reflectflag = 0;
                else
                    m1_ray(t) = m_ray(t-1);
                    m2_ray(t) = m_ray(t-1);
                end
                l1_ray(t) = m1_ray(t).*real(x1_ray(t));
                l2_ray(t) = m2_ray(t).*real(x2_ray(t));
            else
                m1_ray(t) = m_ray(t-5);
                m2_ray(t) = m_ray(t-5);
                l1_ray(t) = m1_ray(t).*real(x1_ray(t-5));
                l2_ray(t) = m2_ray(t).*real(x2_ray(t-5));
                reflectflag = 0;
            end
            
            cgz_1 = (-S2_bilin*l1_ray(t)./(m1_ray(t).^2) - N2_bilin*(l1_ray(t)^2)./(m1_ray(t).^3))/(omega);
            cgz_2 = (-S2_bilin*l2_ray(t)./(m2_ray(t).^2) - N2_bilin*(l2_ray(t)^2)./(m2_ray(t).^3))/(omega);
            cgy_1 = (S2_bilin./m1_ray(t) + N2_bilin*l1_ray(t)./(m1_ray(t)^2))/(omega);
            cgy_2 = (S2_bilin./m2_ray(t) + N2_bilin*l2_ray(t)./(m2_ray(t)^2))/(omega);
            
            % non-hydrostatic group velocities
            %cgz_1 = ((F2_bilin*m_ray(t)*l1_ray(t)^2) + (S2_bilin*l1_ray(t)^3) - ...
            %    (S2_bilin*l1_ray(t)*m_ray(t)^2) - N2_bilin*m_ray(t)*l1_ray(t)^2)/(omega*((m_ray(t)^2 + l1_ray(t)^2)^2));
            %cgz_2 = ((F2_bilin*m_ray(t)*l2_ray(t)^2) + (S2_bilin*l2_ray(t)^3) - ...
            %    (S2_bilin*l2_ray(t)*m_ray(t)^2) - N2_bilin*m_ray(t)*l2_ray(t)^2)/(omega*((m_ray(t)^2 + l2_ray(t)^2)^2));
            %cgy_1 = (-F2_bilin*l1_ray(t)*m_ray(t)^2 + S2_bilin*m_ray(t)^3 - S2_bilin*m_ray(t)*l1_ray(t)^2 + ...
            %    N2_bilin*l1_ray(t)*m_ray(t)^2)/(omega*((m_ray(t)^2 + l1_ray(t)^2)^2));
            %cgy_2 = (-F2_bilin*l2_ray(t)*m_ray(t)^2 + S2_bilin*m_ray(t)^3 - S2_bilin*m_ray(t)*l2_ray(t)^2 + ...
            %    N2_bilin*l2_ray(t)*m_ray(t)^2)/(omega*((m_ray(t)^2 + l2_ray(t)^2)^2));
            
            if ch_ray(t) == 1
                cgz_ray(t) = cgz_1;
                l_ray(t) = l1_ray(t);
                m_ray(t) = m1_ray(t);
                cgy_ray(t) = cgy_1;
                alphaout(t) = x1_ray(t);
            elseif ch_ray(t) == 2
                cgz_ray(t) = cgz_2;
                l_ray(t) = l2_ray(t);
                m_ray(t) = m2_ray(t);
                cgy_ray(t) = cgy_2;
                alphaout(t) = x2_ray(t);
            end
                        
%             % Energy Considerations
%             e_ray(t) = e_ray(1).*sqrt(cgy_ray(1).^2 + cgz_ray(1).^2)./sqrt(cgz_ray(t).^2 + cgy_ray(t).^2);
%             E = @(vs) 0.5.*1025.*(vs.^2 + (-l_ray(t).*vs./m_ray(t)).^2) - e_ray(t);
%             v0 = fzero(E,v_ray(t-1));
            
            %polarization relations
            v_ray(t) = v0;
            w_ray(t) = -l_ray(t).*v_ray(t)./m_ray(t);    % from continuity
            u_ray(t) = v_ray(t).*(F2_bilin + S2_bilin.*l_ray(t)./m_ray(t))./(-1i.*f.*omega);   
            
        end
        
        
    end
    
else
    y_ray(:) = starty.*ones(tsteps,1);
    z_ray(:) = startz.*ones(tsteps,1);
    v_ray(:) = zeros(tsteps,1);
    w_ray(:) = zeros(tsteps,1);
    u_ray(1:tsteps) = zeros(tsteps,1);
    l_ray(1:tsteps) = zeros(tsteps,1);
    m_ray(1:tsteps) = zeros(tsteps,1);
    cgz_ray(1:tsteps) = zeros(tsteps,1);
    cgy_ray(1:tsteps) = zeros(tsteps,1);
end
    



    
end