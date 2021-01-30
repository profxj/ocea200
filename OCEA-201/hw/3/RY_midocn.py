''' %
% This script evaluates the analytical solution for the
% mid-ocean circulation problem of Rhines and Young (1982,
% "A theory for the wind-driven circulation. I. Mid-ocean
% gyres. The domain used corresponds to a square ocean basin
% approximately 5000km wide and is the same as the 2.5 layer
% model used in the companion exercise.
% Follow Pedlosky (1996, Ocean Circulation Theory, Sections
% 3.1-3.7).
%
clear all
close all
'''

import numpy as np

W0 = 1.e-7


def run_midocn(W0=1e-7, ncell=201):
    """
    Run the script

    Parameters
    ----------
    W0 : float
        Choose amplitude of Ekman pumping velocity (m/s).
    ncell : int

    Returns
    -------

    """
    #
    # Set the basin dimensions (m).
    #
    Lx=5.5556e6
    Ly=Lx
    #
    # Set number of points in x and y-directions.
    #
    im=ncell
    jm=ncell
    #
    # Choose the radius of the circular patch of Ekman pumping.
    # (units are in m).
    #
    rekman=1.5e6
    #
    # Set other parameters.
    #
    phi0=40.0
    f0=2.*7.292e-5*np.sin(phi0*2.*4.*np.arctan(1.0)/360.)
    beta=2.0e-11
    #
    # Frictional coupling between layer1 and layer 2.
    #
    A2=1.9e-7
    #
    # Bottom drag coefficient.
    #
    r2=0.0
    #
    # Layer depths. (km??)
    #
    H1=100.0
    H2=100.0
    H=H1+H2
    #
    # Reduced gravities.
    gamma1=7.84e-2
    gamma2=gamma1/2
    #
    # Compute grid-spacing.
    #
    dx=Lx/(im-1)
    dy=Ly/(jm-1)
    #
    # Set the origin at the centre of the domain.
    #
    ix0=(im-1)/2+1
    iy0=(jm-1)/2+1
    x0=(ix0-1)*dx
    y0=(iy0-1)*dy
    #
    # Compute distance of each point from origin.
    #
    r=np.zeros((im,jm))
    x=np.zeros((im,jm))
    y=np.zeros((im,jm))
    for j in range(jm): #1:jm
        for i in range(im): #=1:im
            # Flipping x,y for row,column swap between Python and matlab!
            y[i,j]=(i-1)*dx-x0 + 1
            x[i,j]=(j-1)*dy-y0 + 1
            r[i,j]=np.sqrt(x[i,j]**2 + y[i,j]**2)
    #
    # Compute the forcing (Pedlosky eqn 3.6.1).
    #
    we=np.zeros((im,jm))
    alpha=W0/rekman
    '''
    for j in range(jm):
        for i in range(im):
            if r[i,j]<= rekman
                we[i,j]=-alpha*x[i,j]
    '''
    in_rekman = r < rekman
    we[in_rekman] = -alpha * r[in_rekman]
    #
    # Compute barotropic stream function (Pedlosky eqn 3.6.4).
    #
    psib=np.zeros((im,jm))
    alpha=W0/rekman
    fact=alpha*f0/(2.0*beta*H)

    #for j=1:jm
    # for i=1:im
    #  if r(i,j) <= rekman
    #psib(i,j)=fact*(rekman*rekman-x(i,j)*x(i,j)-y(i,j)*y(i,j));
    psib[in_rekman] =fact*(rekman*rekman-x[in_rekman]**2 -y[in_rekman]**2)
    #
    # Compute geostrophic contours (Pedlosky eqn 3.6.5).
    #
    #q2hat=np.zeros((im,jm))
    fhat=f0*f0*H/(gamma1*H1*H2)
    fact=alpha*f0*fhat/(2.0*beta*H)
    '''
    for j=1:jm
      for i=1:im
       q2hat(i,j)=beta*y(i,j);
       if r(i,j) <= rekman
        q2hat(i,j)=q2hat(i,j)+fact*(rekman*rekman-x(i,j)*x(i,j)-y(i,j)*y(i,j));
       end
      end
    end
    '''
    q2hat = beta * y
    q2hat[in_rekman] = q2hat[in_rekman] + fact * (rekman * rekman - x[in_rekman]**2 -
                                                  y[in_rekman]**2) #y(i, j) * y(i, j));
    #
    # Compute the layer 2 stream function (Pedlosky eqn 3.7.14 and 3.7.15).
    # Psi2 vanishes on the outer-most closed geostrophic contour
    # along which q2hat=q2hatzero=beta*rekman.
    #
    psi2=np.zeros((im,jm))
    fact=A2/(r2*H1/H+A2)
    q2hatzero=beta*rekman
    '''
    for j=1:jm
     for i=1:im
      dq=(q2hat(i,j)-q2hatzero);
      if dq > 0 && y(i,j) < rekman
        psi2(i,j)=fact*dq/fhat;
      end
     end
    end
    '''
    dq = q2hat - q2hatzero
    idx = (dq > 0) & (y < rekman)
    psi2[idx] = fact * dq[idx] / fhat
    #
    # Compute the layer 1 stream function (Pedlosky 3.5.9a).
    #
    psi1=np.zeros((im,jm))
    '''
    for j=1:jm
      for i=1:im
       psi1(i,j)=(H*psib(i,j)-H2*psi2(i,j))/H1;
      end
    end
    '''
    psi1 = (H * psib - H2 * psi2) / H1
    #
    # Compute layer 2 potential vorticty (Pedlosky, eqn 3.7.16).
    #
    G2=f0*f0/(gamma2*H2)
    fact1=((r2*H1/H)-(G2*A2/fhat))/(A2+r2*H1/H)
    fact2=q2hatzero*(1.0+G2/fhat)*A2/(A2+r2*H1/H)
    '''
    for j=1:jm
      for i=1:im
       q2(i,j)=fact1*q2hat(i,j)+fact2;
      end
    end
    '''
    q2 = fact1 * q2hat + fact2

    # Return
    return x, y, r, we, psib, psi1, psi2, q2hat, q2

'''
%
%  Plot the solution.

fsize=6;
%
figure
subplot(3,2,1)
contour(we'); colorbar
title('Ekman pumping velocity, w_E','FontSize',fsize);
subplot(3,2,2)
contour(psib'); colorbar
title('Barotropic stream function, \psi_B','FontSize',fsize);
subplot(3,2,3)
contour(psi1'); colorbar
title('Layer 1 stream function, \psi_1','FontSize',fsize);
subplot(3,2,4)
contour(psi2'); colorbar
title('Layer 2 stream function, \psi_2','FontSize',fsize);
subplot(3,2,5)
contour(q2hat'); colorbar
title('Geostrophic contours','FontSize',fsize);
subplot(3,2,6)
contour(q2'); colorbar
title('Layer 2 potential vorticity, q_2','FontSize',fsize);

print -djpeg -r300 RYfig1.jpg


figure
subplot(3,2,1)
pcolor(we'); shading interp; colorbar
title('Ekman pumping velocity, w_E','FontSize',fsize);
subplot(3,2,2)
pcolor(psib'); shading interp; colorbar
title('Barotropic stream function, \psi_B','FontSize',fsize);
subplot(3,2,3)
pcolor(psi1'); shading interp; colorbar
title('Layer 1 stream function, \psi_1','FontSize',fsize);
subplot(3,2,4)
pcolor(psi2'); shading interp; colorbar
title('Layer 2 stream function, \psi_2','FontSize',fsize);
subplot(3,2,5)
pcolor(q2hat'); shading interp; colorbar
title('Geostrophic contours','FontSize',fsize);
subplot(3,2,6)
pcolor(q2'); shading interp; colorbar
title('Layer 2 potential vorticity, q_2','FontSize',fsize);

print -djpeg -r300 RYfig2.jpg
'''
