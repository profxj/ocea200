""" Port of Andy's Matlab scripts """

import numpy as np

from IPython import embed

#%
#% This script plots the solution for the 2-layer solution of the
#% ventilated thermocline equation of LPS. The script follows
#% Pedlosky (1996, Ocean Circulation Theory), section 4.4.
#%
#% The x-coordinate is longitude in radians, and the y-coordinate
#% is f/f0, starting at the equator.
#%
#%
#%
#%
#%
#%
#%
#%
#%%%%%%%%%%%%%%% DO NOT EDIT THE FILE BELOW HERE %%%%%%%%
#%
def two_layers(theta0=60.,theta2=50., rho1=1025.50,rho2=1026.75,
               Lx=5000., W0=2e-6, H2=400, max_depth=-1200):
    """[summary]

    Args:
        theta0 ([type], optional): [description]. Defaults to 60..
        theta2 ([type], optional): [description]. Defaults to 50..
        rho1 (float, optional): [description]. Defaults to 1025.50.
        rho2 (float, optional): [description]. Defaults to 1026.75.
        Lx ([type], optional): [description]. Defaults to 5000..
        W0 ([type], optional): [description]. Defaults to 2e-6.
        H2 (int, optional): [description]. Defaults to 400.
        max_depth (int, optional): [description]. Defaults to -1200.
    #% Specify the layer densities of the active upper two layers (kg/(m*m*m)).
    #% Northern most extent of the model domain (degrees).
    #% Latitude of the outcrop line for layer 2 (degrees).
    #% Width of the domain (m).
    #% Amplitude of the Ekman pumping velocity (m/s).
    #% Depth of layer 2 along the eastern boundary (m).
    #% NOTE:
    #% Define max plotting depth (m). This parameter controls the maximum value
    #% plotted on the depth axis. You may need to adjust this if in some of your
    #% calculations your layer depths exceed the value of -1200m prescribed here.
    """
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
    theta0=theta0*2*np.pi/360
    theta2=theta2*2*np.pi/360
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
        we[:,j]=-W0*f0*f0*np.sin(np.pi*f[j]/f0)/(f[j]*f[j])
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
            D02[i,j]=D0fact*(1-phi/phie)*np.sin(np.pi*f[j]/f0)
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
    for j in np.where(gdf)[0]:
        for i in range(im): #=1:im
            h[i,j]=np.sqrt((D02[i,j]+H2*H2)/(
                1+gamma1*(1-f[j]/f2)**2/gamma2))
    #
    # Pedlosky eqn 4.4.14a,b.
    #
    h1=np.zeros((im,jm))
    for j in np.where(gdf)[0]:
        for i in range(im): #=1:im
            h1[i,j] = (1-f[j]/f2)*h[i,j]
            h2[i,j] = f[j]*h[i,j]/f2
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
    shadx = np.ones(jm)*phie*eradius/1000
    shady = np.zeros_like(shadx)
    gdj = np.where(gdf)[0]
    phi=np.arange(im)*dphi
    for j in range(jm-1,-1,-1):
        shady[j]=j*dtheta*eradius/1000
        if j in gdj:
            fac=1/(D0fact*np.sin(np.pi*f[j]/f0))
            phi_shadow=phie*(1-fac*gamma1*(1-f[j]/f2)**2*H2*H2/gamma2)
            shadx[j]=phi_shadow*eradius/1000
            #if j == 0:
            #    import pdb; pdb.set_trace()
            gdphi = phi >= phi_shadow
            for i in np.where(gdphi)[0]:
                h[i,j]=H2
                h1[i,j]=np.sqrt(gamma2*D02[i,j]/gamma1)
                h2[i,j]=h[i,j]-h1[i,j]
    import pdb; pdb.set_trace()
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
    D02w=D0fact*np.sin(np.pi*f2/f0)
    Gamma12=gamma1/gamma2
    hw=np.sqrt(D02w+H2*H2)
    pooly=np.arange(jm)*dtheta*eradius/1000
    poolx= np.zeros_like(pooly)

    # Tricky one!
    phi= np.arange(im)*dphi
    for j in np.flip(np.where(gdf)[0]):
        fac=1/(D0fact*np.sin(np.pi*f[j]/f0))
        fac1=Gamma12*(1-f[j]/f2)**2
        phi_pool=phie*(1-fac*(D02w*(1+fac1)+H2*H2*fac1))
        poolx[j] = max(phi_pool*eradius/1000, 0.)
        gdphi = phi <= phi_pool
        for i in np.where(gdphi)[0]:
            h[i,j]=Gamma12*f[j]*hw/(
                f2*(1+Gamma12))+np.sqrt(
                    (D02[i,j]+H2*H2)*(
                        1+Gamma12)-Gamma12*(
                            f[j]*hw/f2)**2)/(1+Gamma12)
            h1[i,j]= h[i,j] - f[j]*hw/f2
            if (i == 10) and (j==10):
                import pdb; pdb.set_trace()
    #
    psi1=np.nan*np.ones((im,jm))
    psi2=np.nan*np.ones((im,jm))

    hp1=h1
    ps=shadx*1000/eradius
    gdf = f <= f2
    for j in np.where(gdf)[0]:
        for i in range(im): #=1:im
            hp1[i,j]=np.nan
            psi2[i,j]=gamma2*h2[i,j]
            psi1[i,j]= gamma1*h1[i,j]+gamma2*(h1[i,j]+h2[i,j])
    import pdb; pdb.set_trace()

    phi=np.arange(im)*dphi
    gdphi = phi < ps
    for j in np.where(gdf)[0]:
        for i in np.where(gdphi)[0]:
            psi2[i,j]=gamma2*(h1[i,j]+h2[i,j])

    # For plotting
    outy=np.ones(jm)*theta2*eradius/1000
    outx=np.arange(im)*dphi*eradius/1000

    ixt = int((xtrans/dx)+1)
    iyt = int(((ytrans*2*np.pi/360)/dtheta)+1)

    return xarr, yarr, shadx, shady, outx, outy, poolx, pooly, psi1, psi2, ixt, iyt, h, hp1

# Command line execution
if __name__ == '__main__':
    two_layers()
    embed(header='225')