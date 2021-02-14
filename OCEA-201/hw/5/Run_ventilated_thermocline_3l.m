%
% This script plots the solution for the 3-layer solution of the
% ventilated thermocline equation of LPS. The script follows
% Pedlosky (1996, Ocean Circulation Theory), section 4.7.
%
% The x-coordinate is longitude in radians, and the y-coordinate
% is f/f0, starting at the equator.
%
clear all
close all
%
% Specify the layer densities of the active upper three layers (kg/(m*m*m)).
%
rho1=????;
rho2=????;
rho3=????;
%
% Northern most extent of model domain (degrees).
%
theta0=??;
%
% Latitude of the outcrop line for layer 3 (degrees).
%
theta3=??;
%
% Latitude of the outcrop line for layer 2 (degrees).
%
theta2=??;
%
% Width of the domain (m).
%
Lx=??????;
%
% Amplitude of the Ekman pumping velocity (m/s).
%
W0=????;
%
% Depth of layer 3 along the eastern boundary (m).
%
H3=???;
%
% NOTE:
% Define max plotting depth (m). This parameter controls the maximum value
% plotted on the depth axis. You may need to adjust this if in some of your
% calculations your layer depths exceed the value of -1200m prescribed here.
%
max_depth=-1200;
%
%%%%%%%%%%%%%%% DO NOT EDIT THE FILE BELOW HERE %%%%%%%%
g=9.81;
rho4=1027.75;
% Layer 1 reduced gravity.
gamma1=(rho2-rho1)*g/rho4;;
% Layer 2 reduced gravity.
gamma2=(rho3-rho2)*g/rho4;;
% Layer 3 reduced gravity.
gamma3=(rho4-rho3)*g/rho4;;
Gamma12=gamma1/gamma2;
Gamma13=gamma1/gamma3;
Gamma23=gamma2/gamma3;
%
% Define grid.
%
im=201;
jm=201;
% Position of y-transect for plotting.
xtrans=Lx/2;
% Position of x-transect for plotting.
ytrans=theta3/2;
% Earth radius.
eradius=6.371e6;
% Angular rotation rate of earth.
Omega=7.292e-5;
theta0=theta0*2*pi/360;
theta3=theta3*2*pi/360;
theta2=theta2*2*pi/360;
f0=2*Omega*sin(theta0);
f3=2*Omega*sin(theta3);
% Latitude grid-spacing.
dtheta=theta0/(jm-1);
j2=int32(theta2/dtheta);
theta2=double(j2)*dtheta;
f2=2*Omega*sin(theta2);
% Longitude grid-spacing.
dx=Lx/(im-1);
dphi=dx/eradius;
phie=(im-1)*dphi;
%
% Coordinate arrays for plotting.
xarr=zeros(im,jm);
yarr=zeros(jm,jm);
for j=1:im
  for i=1:im
   xarr(i,j)=(i-1)*dphi*eradius/1000;
  end
end
for j=1:im
  for i=1:im
   yarr(i,j)=(j-1)*dtheta*eradius/1000;
  end
end
%
% Coriolis parameter.
% 
f=zeros(jm,1);
for j=1:jm
 theta=(j-1)*dtheta;
 f(j)=2*Omega*sin(theta);
end
%
% Ekman pumping - Pedlosky eqn 4.4.25.
%
we=zeros(im,jm);
for j=1:jm
  for i=1:im
    we(i,j)=-W0*f0*f0*sin(pi*f(j)/f0)/(f(j)*f(j));
  end
end
%
% D0^2 from Pedlosky eqn 4.4.26.
% D0^2 from Pedlosky eqn 4.4.26 but NOT using the H3 scaling,
% but instead using the actual factor from 4.4.5 so that H3,
% W0, gamma2, phie and theta0 can be variable parameters.
%
D02=zeros(im,jm);
D0fact=4*eradius*eradius*W0*Omega*sin(theta0)*sin(theta0)*phie/gamma3;
for j=1:jm
  for i=1:im
   phi=(i-1)*dphi;
   D02(i,j)=D0fact*(1-phi/phie)*sin(pi*f(j)/f0);
  end
end
%
% Single layer region f3 <= f <= f0, Pedlosky eqn 4.4.6.
%
h3=zeros(im,jm);
h=zeros(im,jm);
for j=1:jm
  if f(j) >= f3
    for i=1:im
     h3(i,j)=sqrt(D02(i,j)+H3*H3);
     h(i,j)=h3(i,j);
    end
  end
end
%
% Process of subduction, f2 <= f <= f3.
%
% Pedlosky eqn 4.4.18, where h=h2+h3.
%
for j=1:jm
  if f(j) <= f3 & f(j) >= f2
    for i=1:im
     h(i,j)=sqrt((D02(i,j)+H3*H3)/(1+gamma2*(1-f(j)/f3)^2/gamma3));
    end
  end
end
%
% Pedlosky eqn 4.4.14a,b.
%
h2=zeros(im,jm);
for j=1:jm
  if f(j) <= f3 & f(j) >= f2
    for i=1:im
      h2(i,j)=(1-f(j)/f3)*h(i,j);
      h3(i,j)=f(j)*h(i,j)/f3;
    end
  end
end
%
% The layer 3 shadow zone for f2 <= f <= f3.
% The latitude and longitude of the streamline that defines the
% poleward edge of the shadow zone can be computed by equating
% Pedlosky eqn 4.4.26 and 4.4.22.
% Namely:
%  phi=phie*(1-fac*gamma2*(1-f/f3)*(1-f/f3)*H3*H3/gamma3)
% where fac=1/(D0fact*sin(pi*f/f0)).
%
shadx=ones(jm,1)*phie*eradius/1000;
shady=zeros(jm,1);
for j=jm:-1:1
  shady(j)=(j-1)*dtheta*eradius/1000;
  if f(j) < f3 & f(j) >= f2
    for i=1:im
      fac=1/(D0fact*sin(pi*f(j)/f0));
      phi_shadow=phie*(1-fac*gamma2*(1-f(j)/f3)^2*H3*H3/gamma3);
      shadx(j)=phi_shadow*eradius/1000;
      phi=(i-1)*dphi;
      if phi >= phi_shadow
        h(i,j)=H3;
        h2(i,j)=sqrt(gamma3*D02(i,j)/gamma2);
        h3(i,j)=h(i,j)-h2(i,j);
      end
    end
  end
end
%
% The ventilated region solution for layers 1, 2 and 3 (the region
% labelled V in Pedlosky Fig. 4.7.2).
% Use Pelosky eqn 4.7.12 and 4.7.13 to compute h, then use
% eqn 4.7.11a,b,c to compute h1, h2 and h3.
%
h1=zeros(im,jm);
fac2=(1-f2/f3);
for j=1:jm
  if f(j) <= f2
    fac3=(1-f(j)/f3);
    fac4=(1-f(j)/f3-f(j)*fac2*(1+Gamma23*fac3)/(f2*(1+Gamma23*fac2)));
    fac5=fac4*fac4;
    F=1+Gamma23*fac3*fac3+Gamma13*fac5;
    for i=1:im
      h(i,j)=sqrt((D02(i,j)+H3*H3)/F);
      h3(i,j)=f(j)*h(i,j)/f3;
      h2(i,j)=f(j)*h(i,j)*fac2*(1+Gamma23*fac3)/(f2*(1+Gamma23*fac2));
      h1(i,j)=h(i,j)-h2(i,j)-h3(i,j);
    end
  end
end

%
% The layer 3 shadow zone for f <= f2.
% The latitude and longitude of the streamline that defines the
% poleward edge of the shadow zone can be computed by equating
% Pedlosky eqn 4.7.21 and the non-scaled version of 4.4.26.
% Namely:
%  phi=phie(1-fac*(F(f)-1)*h3*h3)
% where fac=1/(D0fact*sin(pi*f/f0)).
%
for j=jm:-1:1
  if f(j) <= f2
    shady(j)=(j-1)*dtheta*eradius/1000;
    fac=1/(D0fact*sin(pi*f(j)/f0));
    fac3=(1-f(j)/f3);
    fac4=(1-f(j)/f3-f(j)*fac2*(1+Gamma23*fac3)/(f2*(1+Gamma23*fac2)));
    fac5=fac4*fac4;
    F=1+Gamma23*fac3*fac3+Gamma13*fac5;
    phi_shadow=phie*(1-fac*(F-1)*H3*H3);
    shadx(j)=phi_shadow*eradius/1000;
    for i=1:im
      phi=(i-1)*dphi;
      if f(j) == f2
        phistar=phi_shadow;
      end
      if phi >= phi_shadow;
        h(i,j)=H3;
%AMM        h2(i,j)=sqrt(gamma3*D02(i,j)/gamma2);
      end
    end
  end
end
%
% The boundary between region M and R in Pedlosky Fig 4.7.2.
% The western boundary of region R is defined by equation 2.40
% if Luyten, Pedlosky and Stommel (1983) by considering a 
% line of constant vorticity in layer 2 which is defined
% by a line of constant (h1+h2) since layer 3 is at rest.
% The constant value of (h1+h2) is given by D02S 
% where D02S=D02(phi^*,theta2).
% The solution in region R is given by Pedlosky eqns 4.7.22,
% 4.7.23 and 4.7.24. While in region M we use 4.7.25 for h2
% with H replaced by H3.
%
RMx=zeros(jm,1);
RMy=zeros(jm,1);
D02S=H3*H3*Gamma23*fac2*fac2;
for j=jm:-1:1
  RMx(j)=NaN;
  RMy(j)=NaN;
  if f(j) <= f2
    RMy(j)=(j-1)*dtheta*eradius/1000;
    for i=1:im
      fac=1/(D0fact*sin(pi*f(j)/f0));
      phi_RM=phie*(1-fac*D02S*(1+Gamma12*(1-f(j)/f2)^2));
      RMx(j)=phi_RM*eradius/1000;
      phi=(i-1)*dphi;
% Solution in region R.
      if phi >= phi_RM
        hhat=sqrt(D02(i,j)/(Gamma23*(1+Gamma12*(1-f(j)/f2)^2)));
        h2(i,j)=f(j)*hhat/f2;
        h1(i,j)=hhat-h2(i,j);
        h3(i,j)=h(i,j)-h1(i,j)-h2(i,j);
      end
% Solution in region M - solution of a quadratic for (h1+h2).
      phis=shadx(j)*1000/eradius;
      if phi >= phis & phi < phi_RM
        a=gamma3*D02(i,j)/gamma1;
        b=f(j)*fac2/(f2*(1+Gamma23*fac2));
        c1=(b*Gamma23-1)^2+gamma2/gamma1;
        c2=2*b*H3*(b*Gamma23-1);
        c3=b*b*H3*H3-a;
        hhat1=(-c2+sqrt(c2*c2-4*c1*c3))/(2*c1);
        hhat2=(-c2-sqrt(c2*c2-4*c1*c3))/(2*c1);
        hhat=hhat1;
%        hhat=hhat2;
        h2(i,j)=b*(H3+gamma2*hhat/gamma3);
        h1(i,j)=hhat-h2(i,j);
        h3(i,j)=h(i,j)-h1(i,j)-h2(i,j);
      end
    end
  end
end

%
%
% The western pool region.
% The latitude and longitude of the streamline that defines the
% eastern edge of the pool region can be found by equating Pedlosky
% eqn 4.6.2 and 4.4.26. It is assumed that the PV is homogenized in the
% pool region which yields Pedlosky eqn 4.6.6 for h and 4.6.5 for h2 in the pool
% in which case h1=h-h2.
% Namely:
% phi=phie*(1-fac*(D02w*(1+gamma1*(1-f/f2)^2/gamma2)
%                 +gamma1*(f-f/f2)^2*H2*H2/(2*gamma2))
% where fac=1/(Dfact*sin(pi*f/f0)), and D02w is the value of D02 evaluated
% at (0,theta2)..
%
poolx=zeros(jm,1);
pooly=zeros(jm,1);
D02w=D0fact*sin(pi*f2/f0);
hw=sqrt(D02w+H3*H3);
for j=jm:-1:1
  pooly(j)=(j-1)*dtheta*eradius/1000;
  if f(j) < f2
    for i=1:im
      fac=1/(D0fact*sin(pi*f(j)/f0));
      fac1=Gamma12*(1-f(j)/f2)^2;
      phi_pool=phie*(1-fac*(D02w*(1+fac1)+H3*H3*fac1));
      poolx(j)=phi_pool*eradius/1000;
      if poolx(j) < 0
        poolx(j)=0;
      end
      phi=(i-1)*dphi;
      if phi <= phi_pool
%AMM       h(i,j)=Gamma12*f(j)*hw/(f2*(1+Gamma12))+...
%AMM              sqrt((D02(i,j)+H3*H3)*(1+Gamma12)-Gamma12*(f(j)*hw/f2)^2)/(1+Gamma12);
%AMM       h1(i,j)=h(i,j)-f(j)*hw/f2;
      end
    end
  end
end

%

%ind1=find(h1==0);
%h1(ind1)=NaN;
%ind2=find(h2==0);
%h2(ind2)=NaN;

hp1=h1;
hp2=h1+h2;
psi1=NaN*ones(im,jm);
psi2=NaN*ones(im,jm);
psi3=NaN*ones(im,jm);
for j=1:jm
  ps=shadx(j)*1000/eradius;
  if f(j) > f3
    for i=1:im
     hp1(i,j)=NaN;
     hp2(i,j)=NaN;
     psi3(i,j)=gamma3*h3(i,j);
    end
  end
  if f(j) <= f3 & f(j) > f2
    for i=1:im
      phi=(i-1)*dphi;
      hp1(i,j)=NaN;
      psi2(i,j)=(gamma2*h2(i,j)+gamma3*(h2(i,j)+h3(i,j)));
      if phi <= ps
        psi3(i,j)=gamma3*(h2(i,j)+h3(i,j));
      end
    end
  end
 if f(j) <= f2
    for i=1:im
      phi=(i-1)*dphi;
      psi1(i,j)=(gamma1*h1(i,j)+gamma2*(h1(i,j)+h2(i,j))+...
                  gamma3*(h1(i,j)+h2(i,j)+h3(i,j)));
      psi2(i,j)=(gamma2*(h1(i,j)+h2(i,j))+gamma3*(h1(i,j)+h2(i,j)+h3(i,j)));
      if phi <= ps
        psi3(i,j)=gamma3*(h1(i,j)+h2(i,j)+h3(i,j));
      end
    end
 end
end


outx2=zeros(im,1);
outy2=ones(jm,1)*theta2*eradius/1000;
for i=1:im
 outx2(i)=(i-1)*dphi*eradius/1000;
end
outx3=zeros(im,1);
outy3=ones(jm,1)*theta3*eradius/1000;
for i=1:im
 outx3(i)=(i-1)*dphi*eradius/1000;
end

ixt=int32(xtrans/dx)+1;
iyt=int32((ytrans*2*pi/360)/dtheta)+1;

figure(1)
subplot(3,2,1)
contour(xarr,yarr,psi3); colorbar
hold on
plot(shadx,shady,'k--')
plot(RMx,RMy,'k--')
plot(outx2,outy2,'k--')
plot(outx3,outy3,'k--')
plot(poolx,pooly,'k--')
xlabel('x (km)')
ylabel('y (km)')
title('a: Layer 3 stream function')

subplot(3,2,3)
contour(xarr,yarr,psi2); colorbar
hold on
plot(shadx,shady,'k--')
plot(RMx,RMy,'k--')
plot(outx2,outy2,'k--')
plot(outx3,outy3,'k--')
plot(poolx,pooly,'k--')
xlabel('x (km)')
ylabel('y (km)')
title('b: Layer 2 stream function')

subplot(3,2,5)
contour(xarr,yarr,psi1); colorbar
hold on
plot(shadx,shady,'k--')
plot(RMx,RMy,'k--')
plot(outx2,outy2,'k--')
plot(outx3,outy3,'k--')
plot(poolx,pooly,'k--')
xlabel('x (km)')
ylabel('y (km)')
title('c: Layer 1 stream function')

subplot(3,2,2)
pcolor(xarr,yarr,psi3); shading interp; colorbar
hold on
plot(shadx,shady,'w--')
plot(RMx,RMy,'w--')
plot(outx2,outy2,'w--')
plot(outx3,outy3,'w--')
plot(poolx,pooly,'w--')
xlabel('x (km)')
ylabel('y (km)')
title('d: Layer 3 stream function')

subplot(3,2,4)
pcolor(xarr,yarr,psi2); shading interp; colorbar
hold on
plot(shadx,shady,'w--')
plot(RMx,RMy,'w--')
plot(outx2,outy2,'w--')
plot(outx3,outy3,'w--')
plot(poolx,pooly,'w--')
xlabel('x (km)')
ylabel('y (km)')
title('e: Layer 2 stream function')

subplot(3,2,6)
pcolor(xarr,yarr,psi1); shading interp; colorbar
hold on
plot(shadx,shady,'w--')
plot(RMx,RMy,'w--')
plot(outx2,outy2,'w--')
plot(outx3,outy3,'w--')
plot(poolx,pooly,'w--')
xlabel('x (km)')
ylabel('y (km)')
title('f: Layer 1 stream function')

print -djpeg -r300 ventilated_3l_fig1_ns.jpg

figure(2)

ixt=int32(xtrans/dx)+1;
iyt=int32((ytrans*2*pi/360)/dtheta)+1;
subplot(3,2,1)
plot(yarr(ixt,:),-h(ixt,:),'b')
hold on
plot(yarr(ixt,:),-hp2(ixt,:),'k')
plot(yarr(ixt,:),-hp1(ixt,:),'r')
axis([0 yarr(ixt,end) max_depth 0]);
xlabel('y (km)')
ylabel('depth (m)')
%legend('h_1+h_2+h_3','h_1+h_2','h_1')
%legend('h_1+h_2+h_3','h_1+h_2','h_1','Location','NorthEastOutside')
title('e: N-S cross-section')

subplot(3,2,3)
plot(xarr(:,iyt),-h(:,iyt),'b')
hold on
plot(xarr(:,iyt),-hp2(:,iyt),'k')
plot(xarr(:,iyt),-hp1(:,iyt),'r')
axis([0 xarr(end,iyt) max_depth 0]);
xlabel('y (km)')
ylabel('depth (m)')
title('f: E-W cross-section')
legend('z_4=h_1+h_2+h_3','z_3=h_1+h_2','z_2=h_1','Location','SouthOutside')

print -djpeg -r300 ventilated_3l_fig2_ns.jpg

figure(3)

subplot(3,2,1)
mesh(xarr,yarr,-h); colorbar
caxis([-1200 0]);
axis([0 xarr(end,iyt) 0 yarr(ixt,end) max_depth 0]);
xlabel('x (km)')
ylabel('y (km)')
zlabel('depth (m)')
title('z_4 surface')


subplot(3,2,3)
mesh(xarr,yarr,-hp2); colorbar
caxis([-1200 0]);
axis([0 xarr(end,iyt) 0 yarr(ixt,end) max_depth 0]);
xlabel('x (km)')
ylabel('y (km)')
zlabel('depth (m)')
title('z_3 surface')

subplot(3,2,5)
mesh(xarr,yarr,-hp1); colorbar
caxis([-1200 0]);
axis([0 xarr(end,iyt) 0 yarr(ixt,end) max_depth 0]);
xlabel('x (km)')
ylabel('y (km)')
zlabel('depth (m)')
title('z_2 surface')

subplot(3,2,2)
mesh(xarr,yarr,-h); colorbar
caxis([-1200 0]);
hold on
mesh(xarr,yarr,-hp2); colorbar
mesh(xarr,yarr,-hp1); colorbar
axis([0 xarr(end,iyt) 0 yarr(ixt,end) max_depth 0]);
xlabel('x (km)')
ylabel('y (km)')
zlabel('depth (m)')
title('z_2, z_3 and z_4 surfaces')

print -djpeg -r300 ventilated_3l_fig3_ns.jpg

% Sverdrup solution test.

svert=gamma3*(D02+H3*H3);
sver=gamma1*h1.*h1+gamma2*(h1+h2).*(h1+h2)+gamma3*h.*h;

figure(4)
subplot(3,2,1)
contour(xarr,yarr,svert); colorbar
xlabel('x (km)')
ylabel('y (km)')
title('Sverdrup solution from D_0^2')
subplot(3,2,3)
contour(xarr,yarr,sver); colorbar
xlabel('x (km)')
ylabel('y (km)')
title('Sverdrup solution from h_1, h_2 and h_3')
subplot(3,2,5)
contour(xarr,yarr,(sver-svert)); colorbar
xlabel('x (km)')
ylabel('y (km)')
title('Sverdrup solution error')

print -djpeg -r300 ventilated_3l_fig4_ns.jpg



