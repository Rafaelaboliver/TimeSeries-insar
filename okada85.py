#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 11:09:30 2017

@author: Phil
"""

#OKADA85 Surface deformation due to a finite rectangular source.
#	[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...
#	   E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN)
#	computes displacements, tilts and strains at the surface of an elastic
#	half-space, due to a dislocation defined by RAKE, SLIP, and OPEN on a 
#	rectangular fault defined by orientation STRIKE and DIP, and size LENGTH and
#	WIDTH. The fault centroid is located (0,0,-DEPTH).
#
#	   E,N    : coordinates of observation points in a geographic referential 
#	            (East,North,Up) relative to fault centroid (meter)
#	   DEPTH  : depth of the fault centroid (DEPTH > 0)
#	   STRIKE : fault trace direction (0 to 360° relative to North), defined so 
#	            that the fault dips to the right side of the trace
#	   DIP    : angle between the fault and a horizontal plane (0 to 90°)
#	   LENGTH : fault length in the STRIKE direction (LENGTH > 0)
#	   WIDTH  : fault width in the DIP direction (WIDTH > 0)
#	   RAKE   : direction the hanging wall moves during rupture, measured relative
#	            to the fault STRIKE (-180 to 180°).
#	   SLIP   : dislocation in RAKE direction (meter)
#	   OPEN   : dislocation in tensile component (meter)
#
#	returns the following variables (same matrix size as E and N):
#	   uN,uE,uZ        : displacements (unit of SLIP and OPEN)
#	   uZE,uZN         : tilts (in rad * FACTOR)
#	   uNN,uNE,uEN,uEE : horizontal strains POSITIVE = COMPRESSION (unit of FACTOR)
#
#	Length unit consistency: E, N, DEPTH, LENGTH, and WIDTH must have the same 
#	unit (e.g. km) which can be different from that of SLIP and OPEN (e.g. m) but
#	with a possible FACTOR on tilt and strain results (in this case, an 
#	amplification of km/m = 1000). To have FACTOR = 1 (tilt in radians and 
#	correct strain unit), use the same length unit for all aforesaid variables.
#
#	Poisson's ratio NU is set to 0.25.
#
#	Formulas and notations from Okada [1985] solution excepted for strain 
#	convention (here positive strain means compression), and for the fault 
#	parameters after Aki & Richards [1980], e.g.:
#	      DIP=90, RAKE=0   : left lateral (senestral) strike slip
#	      DIP=90, RAKE=180 : right lateral (dextral) strike slip
#	      DIP=70, RAKE=90  : reverse fault
#	      DIP=70, RAKE=-90 : normal fault
#
#	Equations are all vectorized excepted for argument DIP which must be
#	a scalar (beacause of a singularity in Okada's equations); all other
#	arguments can be scalar or matrix of the same size.
#
#	Author: Philippe Vernant based on François Beauducel matlab script
#	Created: 2018-01-24
#
#	References:
#	   Aki K., and P. G. Richards, Quantitative seismology, Freemann & Co,
#	      New York, 1980.
#	   Okada Y., Surface deformation due to shear and tensile faults in a
#	      half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.
#
#	Acknowledgments: François Beauducel
#
#   Example :
#       Generate the surface displacement (Ue,Un,Uz) for fault long of 20km 
#       and large of 10km with a N85˚ strike and 30˚ dip with 1m displacement
#       and a rake of 70˚ (left-lateral reverse fault).
#
#       X = np.arange(-50000,50000,5000) # generate a X array from -50km to 50km with a 2km step
#       Y = np.arange(-50000,50000,5000) # generate a Y array from -50km to 50km with a 2km step 
#       X,Y = np.meshgrid(X,Y) # generate a grid from the previous arrays
#       (Ue,Un,Uz)=okada85(X,Y,2500,85,30,20000,10000,90,1,0)
#

# import libraries
from pylab import *
import sys


# Default values for optional input arguments
meshgen   = False    # if E and N must be used to generate a grid
plotmaph  = True  	 # plot map with the fault geometry and the horizontal displacement field
plotmapv  = True	 # plot map with the fault geometry and the vertical displacement field
plotcross = False	 # plot a cross-section profile of the dislocation in the dipping direction
nu = 0.25	         # isotropic Poisson's ratio

def okada85(*args):
    nargin = len(args)
    if nargin < 10 :
        sys.exit("#############################################\n\
                WARNING: Not enough input arguments.\n\
            #############################################\n")
        
    if nargin > 10 :
        sys.exit("#############################################\n\
                WARNING: Too many input arguments.\n\
            #############################################\n")

    # Assigns input arguments
    if meshgen :
        e,n = meshgrid(args[0],args[1])
    else :
        e = args[0]
        n = args[1]
    depth = args[2]
    strike = args[3]*pi/180	# converting STRIKE in radian
    dip = args[4]*pi/180	# converting DIP in radian ('delta' in Okada's equations)
    L = args[5]
    W = args[6]
    rake = args[7]*pi/180	# converting RAKE in radian
    slip = args[8]
    U3 = args[9]

    if dip < 0 or dip > 90 :
        sys.exit("#############################################\n\
                Dip angle must be between 0 and 90.\n\
            #############################################\n")
        
    if round((W/2*sin(dip)),10)  > depth :
        sys.exit("#############################################\n\
                The dislocation top edge is above the ground surface.\n\
            #############################################\n")
#    else :
#        print('the dislocation top edge is %0.3fm below the surface'%(depth-W/2*sin(dip)))
        
    # Defines dislocation in the fault plane system
    U1 = cos(rake)*slip
    U2 = sin(rake)*slip
    
    # Converts fault coordinates (E,N,DEPTH) relative to centroid
    # into Okada's reference system (X,Y,D)
    d = depth + sin(dip)*W/2	# d is fault's top edge
    ec = e + cos(strike)*cos(dip)*W/2
    nc = n - sin(strike)*cos(dip)*W/2
    x = cos(strike)*nc + sin(strike)*ec + L/2
    y = sin(strike)*nc - cos(strike)*ec + cos(dip)*W
    
    # Variable substitution (independent from xi and eta)
    p = y*cos(dip) + d*sin(dip)
    q = y*sin(dip) - d*cos(dip)    

    # Displacements
    #U1 () strike-slip component
    #U2 () dip-slip component
    #U3 () tensile fault component
    

    ux = -U1/(2*pi) * ( ux_ss(x,p,q,dip,nu) - ux_ss(x,p-W,q,dip,nu) \
              - ux_ss(x-L,p,q,dip,nu) + ux_ss(x-L,p-W,q,dip,nu))    \
         -U2/(2*pi) * (ux_ds(x,p,q,dip,nu) - ux_ds(x,p-W,q,dip,nu)  \
              - ux_ds(x-L,p,q,dip,nu) + ux_ds(x-L,p-W,q,dip,nu))    \
         +U3/(2*pi) * (ux_tf(x,p,q,dip,nu) - ux_tf(x,p-W,q,dip,nu)  \
              - ux_tf(x-L,p,q,dip,nu) + ux_tf(x-L,p-W,q,dip,nu))            
             
    uy = -U1/(2*pi) * ( uy_ss(x,p,q,dip,nu) - uy_ss(x,p-W,q,dip,nu) \
              - uy_ss(x-L,p,q,dip,nu) + uy_ss(x-L,p-W,q,dip,nu))    \
         -U2/(2*pi) * (uy_ds(x,p,q,dip,nu) - uy_ds(x,p-W,q,dip,nu)  \
              - uy_ds(x-L,p,q,dip,nu) + uy_ds(x-L,p-W,q,dip,nu))    \
         +U3/(2*pi) * (uy_tf(x,p,q,dip,nu) - uy_tf(x,p-W,q,dip,nu)  \
              - uy_tf(x-L,p,q,dip,nu) + uy_tf(x-L,p-W,q,dip,nu))            

    uz = -U1/(2*pi) * ( uz_ss(x,p,q,dip,nu) - uz_ss(x,p-W,q,dip,nu) \
              - uz_ss(x-L,p,q,dip,nu) + uz_ss(x-L,p-W,q,dip,nu))    \
         -U2/(2*pi) * (uz_ds(x,p,q,dip,nu) - uz_ds(x,p-W,q,dip,nu)  \
              - uz_ds(x-L,p,q,dip,nu) + uz_ds(x-L,p-W,q,dip,nu))    \
         +U3/(2*pi) * (uz_tf(x,p,q,dip,nu) - uz_tf(x,p-W,q,dip,nu)  \
              - uz_tf(x-L,p,q,dip,nu) + uz_tf(x-L,p-W,q,dip,nu))            

    # Rotation from Okada's axes to geographic
    ue = sin(strike) * ux - cos(strike) * uy
    un = cos(strike) * ux + sin(strike) * uy

# =================================================================
# plots
    
    alpha = pi/2 - strike
    x_fault = L/2*cos(alpha)*array([-1,1,1,-1]) + sin(alpha)*cos(dip)*W/2*array([-1,-1,1,1])
    y_fault = L/2*sin(alpha)*array([-1,1,1,-1]) + cos(alpha)*cos(dip)*W/2*array([1,1,-1,-1])
    z_fault = -d + sin(dip)*W*array([1,1,0,0])
    ddx = U1*cos(alpha) - U2*sin(alpha)*cos(dip) + U3*sin(alpha)*sin(dip)
    ddy = U1*sin(alpha) + U2*cos(alpha)*cos(dip) - U3*cos(alpha)*sin(dip)
    ddz = U2*sin(dip) + U3*cos(dip)

    if plotmaph:
        # 2D plot
        figure(201,figsize=(10,10))
        clf()
        contourf(e,n,np.sqrt(ue**2+un**2))  
        cbar=colorbar(shrink=0.5)
        cbar.ax.set_ylabel('horizontal displacement')
        plot([0],[0],'rx',label="centroid")
        quiver([0],[0],ddx,ddy,color='r',units='width')
        plot(append(x_fault,x_fault[0]),append(y_fault,y_fault[0]),label="fault plane")
        plot(x_fault[0:2],y_fault[0:2],'r-',label="top fault plane edge")
        plot(e,n,'k.')
        quiver(e,n,ue,un,units='width')
        axis('scaled')
        xlabel('East')
        ylabel('North')
        legend(loc='upper right')
        title("Horizontal displacement")
        savefig('horizontal.png',dpi=300)

    if plotmapv:
        # 2D plot
        figure(202,figsize=(10,10))
        clf()
        contourf(e,n,uz) 
        cbar=colorbar(shrink=0.5)
        cbar.ax.set_ylabel('vertical displacement')
        scatter(e, n, s=1)
        plot([0],[0],'rx',label="centroid")
        plot(append(x_fault,x_fault[0]),append(y_fault,y_fault[0]),label="fault plane")
        plot(x_fault[0:2],y_fault[0:2],'r-',label="top fault plane edge")
        quiver([0],[0],ddx,ddy,color='r',units='width')
        axis('scaled')
        xlabel('East')
        ylabel('North')
        legend(loc='upper right')
        title("Vertical displacement")
        savefig('vertical.png',dpi=300)

    if plotcross:
        # 2D cross section
        figure(203)
        clf()
        xp = [-W/2 * cos(dip) , W/2 * cos(dip)] 
        yp = [-depth + W/2 * sin(dip) , -depth - W/2 * sin(dip)] 
        plot(xp,yp,'b-',label="fault plane")
        plot(0,-depth,'ro',label="centroid")
        legend(loc='upper right')
#        axis('scaled')
        xlabel('Distance normal to strike')
        ylabel('Depth')
        savefig('cross_section.png',dpi=300)

#    np.savetxt('okada.txt',np.c_[reshape(e,(size(e),1)),reshape(n,(size(e)),1),\
#            reshape(ue,(size(e),1)),reshape(un,(size(e),1)),reshape(uz,(size(e),1))])  # ecriture du fichier

    return ue,un,uz

# =================================================================
# Displacement subfunctions

# strike-slip displacement subfunctions [equation (25) p. 1144]
# -----------------------------------------------------------------
def ux_ss (xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    u = xi*q/(R*(R + eta)) + I1(xi,eta,q,dip,nu,R)*sin(dip)
    u = where(q!=0,u+arctan(xi*eta/(q*R)),u)

    return u

# -----------------------------------------------------------------
def uy_ss(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    u = (eta*cos(dip) + q*sin(dip))*q/(R*(R + eta)) \
        	+ q*cos(dip)/(R + eta) + I2(eta,q,dip,nu,R)*sin(dip)

    return u

# -----------------------------------------------------------------
def uz_ss(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(dip) - q*cos(dip)
    u = (eta*sin(dip) - q*cos(dip))*q/(R*(R + eta)) \
        + q*sin(dip)/(R + eta) + I4(db,eta,q,dip,nu,R)*sin(dip)

    return u

# dip-slip displacement subfunctions [equation (26) p. 1144]
# -----------------------------------------------------------------
def ux_ds(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q **2)
    u = q/R - I3(eta,q,dip,nu,R)*sin(dip)*cos(dip)

    return u

# -----------------------------------------------------------------
def uy_ds(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    u = (eta*cos(dip) + q*sin(dip))*q/(R*(R + xi)) \
        - I1(xi,eta,q,dip,nu,R)*sin(dip)*cos(dip)
    u = where(q!=0,u+cos(dip)*arctan(xi*eta/(q*R)),u)

    return u

# -----------------------------------------------------------------
def uz_ds(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(dip) - q*cos(dip)
    u = db*q/(R*(R + xi)) \
        - I5(xi,eta,q,dip,nu,R,db)*sin(dip)*cos(dip)
    u = where(q!=0,u+sin(dip)*arctan(xi*eta/(q*R)),u)

    return u

# tensile fault displacement subfunctions [equation (27) p. 1144]
# -----------------------------------------------------------------
def ux_tf(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    u = q**2 /(R*(R + eta)) \
        - I3(eta,q,dip,nu,R)*sin(dip)**2

    return u

# -----------------------------------------------------------------
def uy_tf(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    u = -(eta*sin(dip) - q*cos(dip))*q/(R*(R + xi)) \
        - sin(dip)*xi*q/(R*(R + eta)) \
        - I1(xi,eta,q,dip,nu,R)*sin(dip)**2
    u = where(q!=0,u+sin(dip)*arctan(xi*eta/(q*R)),u)

    return u

# -----------------------------------------------------------------
def uz_tf(xi,eta,q,dip,nu) :
    R = sqrt(xi**2 + eta**2 + q**2)
    db = eta*sin(dip) - q*cos(dip)
    u = (eta*cos(dip) + q*sin(dip))*q/(R*(R + xi)) \
        + cos(dip)*xi*q/(R*(R + eta)) \
        - I5(xi,eta,q,dip,nu,R,db)*sin(dip)**2
    u = where(q!=0,u-cos(dip)*arctan(xi*eta/(q*R)), u)

    return u

# I... displacement subfunctions [equations (28) (29) p. 1144-1145]
# -----------------------------------------------------------------
def I1(xi,eta,q,dip,nu,R) :
    db = eta*sin(dip) - q*cos(dip)
    if cos(dip) > spacing(1) :
        I = (1 - 2*nu) * (-xi/(cos(dip)*(R + db))) \
            - sin(dip)/cos(dip)*I5(xi,eta,q,dip,nu,R,db)
    else :
        I = -(1 - 2*nu)/2 * xi*q/(R + db)**2

    return I

# -----------------------------------------------------------------
def I2(eta,q,dip,nu,R) :
    I = (1 - 2*nu) * (-log(R + eta)) - I3(eta,q,dip,nu,R)

    return I

# -----------------------------------------------------------------
def I3(eta,q,dip,nu,R) :
    yb = eta*cos(dip) + q*sin(dip)
    db = eta*sin(dip) - q*cos(dip)
    if cos(dip) > spacing(1) :
        I = (1 - 2*nu) * (yb/(cos(dip)*(R + db)) - log(R + eta)) \
            + sin(dip)/cos(dip) * I4(db,eta,q,dip,nu,R)
    else :
        I = (1 - 2*nu)/2 * (eta/(R + db) + yb*q/(R + db)**2 - log(R + eta))

    return I

# -----------------------------------------------------------------
def I4(db,eta,q,dip,nu,R) : 
    if cos(dip) > spacing(1) : 
        I = (1 - 2*nu) * 1/cos(dip) * (log(R + db) - sin(dip)*log(R + eta))
    else :
        I = -(1 - 2*nu) * q/(R + db)

    return I

# -----------------------------------------------------------------
def I5(xi,eta,q,dip,nu,R,db) :
    X = sqrt(xi**2 + q**2)
    if cos(dip) > spacing(1) :
        I = (1 - 2*nu) * 2./cos(dip) * arctan((eta*(X + q*cos(dip)) \
            + X*(R + X)*sin(dip)) /(xi*(R + X)*cos(dip)))
        I = where(xi==0,0,I)
    else:
        I = -(1 - 2*nu) * xi*sin(dip)/(R + db)

    return I


