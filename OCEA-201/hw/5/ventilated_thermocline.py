""" Port of Andy's Matlab scripts """

import numpy as np

#%
#% This script plots the solution for the 2-layer solution of the
#% ventilated thermocline equation of LPS. The script follows
#% Pedlosky (1996, Ocean Circulation Theory), section 4.4.
#%
#% The x-coordinate is longitude in radians, and the y-coordinate
#% is f/f0, starting at the equator.
#%
#% Specify the layer densities of the active upper two layers (kg/(m*m*m)).
rho1=1#????;
rho2=1#????;
#%
#% Northern most extent of the model domain (degrees).
#%
theta0=45#??;
#%
#% Latitude of the outcrop line for layer 2 (degrees).
#%
theta2=45#??;
#%
#% Width of the domain (m).
#%
Lx=1000#???????;
#%
#% Amplitude of the Ekman pumping velocity (m/s).
#%
W0=0.1#????;
#%
#% Depth of layer 2 along the eastern boundary (m).
#%
H2=100#???;
#%
#% NOTE:
#% Define max plotting depth (m). This parameter controls the maximum value
#% plotted on the depth axis. You may need to adjust this if in some of your
#% calculations your layer depths exceed the value of -1200m prescribed here.
#%
max_depth=-1200
#%
#%%%%%%%%%%%%%%% DO NOT EDIT THE FILE BELOW HERE %%%%%%%%
#%
def two_layers():
    g=9.81
    rho3=1027.50
    #%
    # Layer 1 reduced gravity.
    gamma1=(rho2-rho1)*g/rho3
    # Layer 2 reduced gravity.
    gamma2=(rho3-rho2)*g/rho3
    #
    # Define grid.
    #
    im=201
    jm=201
    # Position of y-transect for plotting.
    xtrans=Lx/2
    # Position of x-transect for plotting.
    ytrans=theta0/2
    # Earth radius.
    eradius=6.371e6

    # Angular rotation rate of earth.
    Omega=7.292e-5
    theta0=theta0*2*pi/360
    theta2=theta2*2*pi/360
    f0=2*Omega*np.sin(theta0)
    f2=2*Omega*np.sin(theta2)
    # Latitude grid-spacing.
    dtheta=theta0/(jm-1)
    # Longitude grid-spacing.
    dx=Lx/(im-1)
    dphi=dx/eradius
    phie=(im-1)*dphi
    #
    # Coordinate arrays for plotting.
    xarr=np.zeros((im,jm))
    yarr=np.zeros((jm,jm))

    for i in range(im): #=1:im
        xarr[i,:]=(i-1)*dphi*eradius/1000
    for j in range(jm): #1:jm
        yarr[i,j]=(j-1)*dtheta*eradius/1000
    #
    # Coriolis parameter.
    # 
    #for j=1:jm
    theta= np.arange(jm)*dtheta
    f=2*Omega*np.sin(theta)

    #
    # Ekman pumping - Pedlosky eqn 4.4.25.
    #
    we=np.zeros((im,jm))
    for j in range(jm): #1:jm
        we[:,j]=-W0*f0*f0*np.sin(np.pi*f[j]/f0)/(f[j]*f[j]);
    #
    # D0^2 from Pedlosky eqn 4.4.26 but NOT using the H2 scaling,
    # but instead using the actual factor from 4.4.5 so that H2,
    # W0, gamma2, phie and theta0 can be variable parameters.
    #
    D02=np.zeros((im,jm))
    D0fact=4*eradius*eradius*W0*Omega*np.sin(
        theta0)*np.sin(theta0)*phie/gamma2

    for j in range(jm): #1:jm
        for i in range(im): #=1:im
            phi=i*dphi
            D02[i,j]=D0fact*(1-phi/phie)*np.sin(pi*f[j]/f0)
    #
    # Single layer region f0 <= f <= f2, Pedlosky eqn 4.4.6.
    #
    #   h2(i,j)=sqrt(D02(i,j)+H2*H2);
    #   h(i,j)=h2(i,j);
    h2 = np.sqrt(D02+H2*H2)
    h = h2.copy()
    #
    # Process of subduction, f2 < f <= 0..
    #
    # Pedlosky eqn 4.4.18, where h=h1+h2.
    #
    #for j=1:jm
    #if f(j) <= f2
    #    for i=1:im
    #    h(i,j)=sqrt((D02(i,j)+H2*H2)/(1+gamma1*(1-f(j)/f2)^2/gamma2));
    #    end
    
    gdf = f <= f2
    for i in range(im): #=1:im
        h[i,gdf]=np.sqrt((D02[i,gdf]+H2*H2)/(
            1+gamma1*(1-f[gdf]/f2)**2/gamma2))

    #
    # Pedlosky eqn 4.4.14a,b.
    #
    h1=np.zeros((im,jm))
    for i in range(im): #=1:im
        h1[i,gdf] = (1-f[gdf]/f2)*h[i,gdf]
        h2[i,gdf] = f[gdf]*h[i,gdf]/f2
    #
    # The shadow zone.
    # The latitude and longitude of the streamline that defines the
    # poleward edge of the shadow zone can be computed by equating
    # Pedlosky eqn 4.4.26 and 4.4.22.
    # Namely:
    #  phi=phie*(1-fac*gamma1*(1-f/f2)*(1-f/f2)*H2*H2/gamma2)
    # where fac=1/(D0fact*sin(pi*f/f0)).
    #
    #shadx=ones(jm,1)*phie*eradius/1000;
    #shady=zeros(jm,1);
    #for j=jm:-1:1
    shady = np.arange(jm)*dtheta*eradius/1000
    shadx = np.zeros_like(shady)
    j = gdf
    fac=1/(D0fact*np.sin(np.pi*f[j]/f0))
    phi_shadow=phie*(1-fac*gamma1*(1-f[j]/f2)**2*H2*H2/gamma2)
    shadx[j]=phi_shadow*eradius/1000
    phi=np.arange(im)*dphi
    gdphi = phi >= phi_shadow
    for j in np.where(gdf)[0]:
        for i in np.where(gdphi)[0]:
            h[i,j]=H2
            h1[i,j]=np.sqrt(gamma2*D02[i,j]/gamma1)
            h2[i,j]=h[i,j]-h1[i,j]
    #
    # The western pool region.
    # The latitude and longitude of the streamline that defines the
    # eastern edge of the pool region can be found by equating Pedlosky
    # eqn 4.6.2 and 4.4.26. It is assumed that the PV is homogenized in the
    # pool region which yields Pedlosky eqn 4.6.6 for h and 4.6.5 for h2 in the pool
    # in which case h1=h-h2.
    # Namely:
    # phi=phie*(1-fac*(D02w*(1+gamma1*(1-f/f2)^2/gamma2)/(2*H2^2)
    #                 +gamma1*(f-f/f2)^2/(2*gamma2))
    # where fac=1/(D0fact*sin(pi*f/f0)), and D02w is the value of D02 evaluated
    # at (0,theta2)..
    #
    poolx=zeros(jm,1);
    D02w=D0fact*np.sin(pi*f2/f0)
    Gamma12=gamma1/gamma2
    hw=np.sqrt(D02w+H2*H2)
    pooly=np.arange(jm)*dtheta*eradius/1000

    # Tricky one!
    fac=1/(D0fact*np.sin(pi*f[gdf]/f0))
    fac1=Gamma12*(1-f[gdf]/f2)^2
    phi_pool=phie*(1-fac*(D02w*(1+fac1)+H2*H2*fac1))
    poolx= np.maximum(phi_pool*eradius/1000, 0.)
    phi= np.arange(im)*dphi
    gdphi = phi <= phi_pool
    for j in np.where(gdf)[0]:
        for i in np.where(gdphi)[0]:
            h[i,j]=Gamma12*f[j]*hw/(
                f2*(1+Gamma12))+np.sqrt(
                    (D02[i,j]+H2*H2)*(
                        1+Gamma12)-Gamma12*(
                            f[j]*hw/f2)**2)/(1+Gamma12)
            h1[i,j]=h[i,j]-f[j]*hw/f2
    #
    psi2=np.nan*np.ones((im,jm))

    hp1=h1
    ps=shadx*1000/eradius
    for i in range(im): #=1:im
        hp1[i,gdf]=np.nan
        psi2[i,gdf]=gamma2*h2[i,gdf]

    phi=np.arange(im)*dphi
    psi1= gamma1*h1+gamma2*(h1+h2)
    gdphi = phi < ps
    for j in np.where(gdf)[0]:
        for i in np.where(gdphi)[0]:
            psi2[i,j]=gamma2*(h1[i,j]+h2[i,j])

    # For plotting
    outy=np.ones(jm)*theta2*eradius/1000
    outx=(np.arange(im)*dphi*eradius/1000

    '''
    figure(1)
    subplot(3,2,1)
    contour(xarr,yarr,psi2); colorbar
    hold on
    plot(shadx,shady,'k--')
    plot(outx,outy,'k--')
    plot(poolx,pooly,'k--')
    xlabel('x (km)')
    ylabel('y (km)')
    title('a: Layer 2 stream function')
    jlab=50;
    %text(jlab*dphi*eradius/1000,0.7*shady(jlab),'Shadow Zone')

    subplot(3,2,3)
    contour(xarr,yarr,psi1); colorbar
    hold on
    plot(shadx,shady,'k--')
    plot(outx,outy,'k--')
    plot(poolx,pooly,'k--')
    xlabel('x (km)')
    ylabel('y (km)')
    title('b: Layer 1 stream function')

    subplot(3,2,2)
    pcolor(xarr,yarr,psi2); shading interp; colorbar
    hold on
    plot(shadx,shady,'w--')
    plot(outx,outy,'w--')
    plot(poolx,pooly,'w--')
    xlabel('x (km)')
    ylabel('y (km)')
    title('c: Layer 2 stream function')
    %text(jlab*dphi*eradius/1000,0.7*shady(jlab),'Shadow Zone')

    subplot(3,2,4)
    pcolor(xarr,yarr,psi1); shading interp; colorbar
    hold on
    plot(shadx,shady,'w--')
    plot(outx,outy,'w--')
    plot(poolx,pooly,'w--')
    xlabel('x (km)')
    ylabel('y (km)')
    title('d: Layer 1 stream function')

    ixt=int32(xtrans/dx)+1;
    iyt=int32((ytrans*2*pi/360)/dtheta)+1;
    subplot(3,2,5)
    plot(yarr(ixt,:),-h(ixt,:),'b')
    hold on
    plot(yarr(ixt,:),-hp1(ixt,:),'r')
    axis([0 yarr(ixt,end) max_depth 0]);
    xlabel('y (km)')
    ylabel('depth (m)')
    title('e: N-S cross-section')
    legend('z_3=h_1+h_2','z_2=h_1')

    subplot(3,2,6)
    plot(xarr(:,iyt),-h(:,iyt),'b')
    hold on
    plot(xarr(:,iyt),-hp1(:,iyt),'r')
    axis([0 xarr(end,iyt) max_depth 0]);
    xlabel('y (km)')
    ylabel('depth (m)')
    title('f: E-W cross-section')
    %legend('z_3=h_1+h_2','z_2=h_1','Location','SouthOutside')

    print -djpeg -r300 ventilated_2l_fig1_ns.jpg

    figure(2)

    subplot(3,2,1)
    mesh(xarr,yarr,-h); colorbar
    caxis([-1200 0]);
    axis([0 xarr(end,iyt) 0 yarr(ixt,end) max_depth 0]);
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('depth (m)')
    title('z_3 surface')


    subplot(3,2,3)
    mesh(xarr,yarr,-hp1); colorbar
    caxis([-1200 0]);
    axis([0 xarr(end,iyt) 0 yarr(ixt,end) max_depth 0]);
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('depth (m)')
    title('z_2 surface')

    subplot(3,2,5)
    mesh(xarr,yarr,-h); colorbar
    caxis([-1200 0]);
    hold on
    mesh(xarr,yarr,-hp1); colorbar
    axis([0 xarr(end,iyt) 0 yarr(ixt,end) max_depth 0]);
    xlabel('x (km)')
    ylabel('y (km)')
    zlabel('depth (m)')
    title('z_2 and z_3 surfaces')

    print -djpeg -r300 ventilated_2l_fig2_ns.jpg
    '''

