# Arya p198, eq 9.3

# c = Q/(2 * pi * ubar * sigy * sigz) 
#   * exp( -y^2/(2 * sigy^2) ) 
#   * ( 
#       exp( -(z-H)^2/(2 * sigz^2) ) 
#     + exp( -(z+H)^2/(2 * sigz^2) ) 
#     )

import numpy as np
import warnings

# for sigy
c = np.array([24.1667, 18.333, 12.5  , 8.333  , 6.25   , 4.1667  ])
d = np.array([ 2.5334, 1.8096, 1.0857, 0.72382, 0.54287, 0.36191 ])

nzbrur = np.array([8,3,1,6,9,10])
xzbrur = np.array([
    0.10,   0.15,   0.20,   0.25,   0.30,   0.40,   0.50,   3.11,   np.inf, np.inf,
    0.20,   0.40,   np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
    np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 
    0.30,   1.00,   3.00,   10.0,   30.0,   np.inf, np.inf, np.inf, np.inf, np.inf, 
    0.10,   0.30,   1.00,   2.00,   4.00,   10.0,   20.0,   40.0,   np.inf, np.inf, 
    0.20,   0.70,   1.00,   2.00,   3.00,   7.00,   15.0,   30.0,   60.0,   np.inf,
    ]).reshape((6, 10))
azrur = np.array([
    122.80, 158.08, 170.22, 179.52, 217.41, 258.89, 346.75, 453.85, np.nan, np.nan,
    90.673, 98.483, 109.3,  np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
    61.141, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
    34.459, 32.093, 32.093, 33.504, 36.65 , 44.053, np.nan, np.nan, np.nan, np.nan, 
    24.26 , 23.331, 21.628, 21.628, 22.534, 24.703, 26.97 , 35.42 , 47.618, np.nan,
    15.209, 14.457, 13.953, 13.953, 14.823, 16.187, 17.836, 22.651, 27.074, 34.219,
    ]).reshape((6, 10))
bzrur = np.array([
    0.9447 , 1.0542 , 1.0932 , 1.1262 , 1.2644 , 1.4094 , 1.7283 , 2.1166 , np.nan , np.nan ,
    0.93198, 0.98332, 1.0971 , np.nan , np.nan , np.nan , np.nan , np.nan , np.nan , np.nan , 
    0.91465, np.nan , np.nan , np.nan , np.nan , np.nan , np.nan , np.nan , np.nan , np.nan , 
    0.86974, 0.81066, 0.64403, 0.60486, 0.56589, 0.51179, np.nan , np.nan , np.nan , np.nan , 
    0.8366 , 0.81956, 0.75660, 0.63077, 0.57154, 0.50527, 0.46713, 0.37615, 0.29592, np.nan ,
    0.81558, 0.78407, 0.68465, 0.63227, 0.54503, 0.46490, 0.41507, 0.32681, 0.27436, 0.21716,
    ]).reshape((6, 10))
szbrur = azrur * (xzbrur ** bzrur)
for i,j in enumerate(nzbrur):
    szbrur[i, j-1 ] = np.inf


def pasquil_gifford_z(istab, dxkm, irus=0, mrough=0, sigz0=0.001):

    istab0 = istab - 1  # 0-5 instead of 1-6
    if irus != 0: raise NotImplemented('non rural')

    zerosig = .001
    # blockdata.f:
    # data nzbrur/8,3,1,6,9,10/

    # csigma.puf:
    #              azrur(10,6),bzrur(10,6),
    # blockdata.f:
    #   data xzbrur/0.10, 0.15, 0.20, 0.25, 0.30,
    #  1            0.40, 0.50, 3.11, 2*1.e20,
    #  2            0.20, 0.40, 8*1.e20,
    #  3            10*1.e20,
    #  4            0.30, 1.00, 3.00, 10.0, 30.0, 5*1.e20,
    #  5            0.10, 0.30, 1.00, 2.00, 4.00,
    #  5            10.0, 20.0, 40.0, 2*1.e20,
    #  6            0.20, 0.70, 1.00, 2.00, 3.00,
    #  6            7.00, 15.0, 30.0, 60.0, 1.e20/
    #  data azrur/122.80, 158.08, 170.22, 179.52, 217.41,
    # 1           258.89, 346.75, 3*453.85,
    # 2           90.673, 98.483, 8*109.3,
    # 3           10*61.141,
    # 4           34.459, 32.093, 32.093, 33.504, 36.65 ,5*44.053,
    # 5           24.26 , 23.331, 21.628, 21.628, 22.534,
    # 5           24.703, 26.97 , 35.42 , 2*47.618,
    # 6           15.209, 14.457, 13.953, 13.953, 14.823,
    # 6           16.187, 17.836, 22.651, 27.074, 34.219/
    #  data bzrur/0.9447 , 1.0542 , 1.0932 , 1.1262 , 1.2644 ,
    # 1           1.4094 , 1.7283,  3*2.1166,
    # 2           0.93198, 0.98332, 8*1.0971,
    # 3           10*0.91465,
    # 4           0.86974, 0.81066, 0.64403, 0.60486, 0.56589, 5*0.51179,
    # 5           0.8366 , 0.81956, 0.75660, 0.63077, 0.57154,
    # 5           0.50527, 0.46713, 0.37615, 2*0.29592,
    # 6           0.81558, 0.78407, 0.68465, 0.63227, 0.54503,
    # 6           0.46490, 0.41507, 0.32681, 0.27436, 0.21716/

    # sigset.f:
    # c     Additive term to adjust B-exponent of sigma-z
    #       bz0trm(i,j,md)=0.0777+.0215*alog(z0(i,j,md))
    # c     Factor to adjust A-coefficient of sigma-z
    #       az0fac(i,j,md)=(1000.**bz0trm(i,j,md))*
    #&                      1.585*z0(i,j,md)**0.1301
    # c     PG and MESOPUFF sigma-y factor
    #       yz0fac(i,j,md)=(z0(i,j,md)/.03)**0.2
    #       ymaxfac=MAX(ymaxfac,yz0fac(i,j,md))
    # csigma.f:
    # c        IRUS - integer   - Rural cell indicator (rural=0 ; urban=1)


    ############################
    # sigtz.f
    ############################
    # c  Compute rural PGT sigma z (in m) using the ISC formulation.
    #    if(idoptzs.eq.3 .and. irus.eq.0) then
 
    # c     Define the number of sigma Z curve segments, NSEG.
    #       nseg = nzbrur(istab)

    # c     Set the roughness-adjusted A', B' values for this location
    # c     and stability class
    #       aprm(:)=azrur(:,istab)
    #       bprm(:)=bzrur(:,istab)
    #       if(bz0trm(ixcell,iycell,mdcell).GT.-999.0) then
    # c        Apply roughness adjustment (B=-1000. if MROUGH=0)
    #          do k=1,nseg
    #             aprm(k)=azrur(k,istab)*az0fac(ixcell,iycell,mdcell)
    #             bprm(k)=bzrur(k,istab)-bz0trm(ixcell,iycell,mdcell)
    #          enddo
    #       endif
    # c     Set sigma-z values at segment boundaries
    #       do k=1,nseg-1
    #          szbrur(k)=aprm(k)*(xzbrur(k,istab))**bprm(k)
    #       enddo
    #       szbrur(nseg)=1.0e20
 
    #       xpskm = zero
    # c     If initial sigma, sigz0, solve for pseudo distance, xpskm
    #       if(sigz0.gt.zerosig) then
 
    # c        Find the correct region, denoted by the segment JSEG.
    #          jseg = 1
    #          if(sigz0.le.szbrur(jseg)) go to 20
    #          jseg = jseg + 1
    #          if(jseg.lt.nseg) go to 15
 
    #          xpskm = ( sigz0/aprm(jseg) )**(1./bprm(jseg))
    #       endif
 
    # c     Combine the actual and pseudo distances and compute sigz.
    #       xkm = xpskm + dxkm
 
    # c     Compute sigma
    #       if(dxkm.GT.zero) then
    # c        Find the correct region, denoted by the segment JSEG.
    #          jseg = 1
    #          if(xkm.le.xzbrur(jseg,istab)) go to 30
    #          jseg = jseg + 1
    #          if(jseg.lt.nseg) go to 25
 
    #          sigz = aprm(jseg) * xkm**bprm(jseg)
    #       else
    #          sigz = sigz0
    #       endif
 
    # c     Assign virtuals
    #       virdz = xkm
    #       virtz = xkm *1000. / uavgs
 
    #    endif


    if mrough != 0: 
        # TODO aprm, bprm get adjusted
        raise NotImplemented('mrough')


    assert dxkm.shape == istab0.shape
    sigz = np.empty(dxkm.shape, dtype=float)

    for i, (dx, st) in enumerate(zip(np.array(dxkm).flat, np.array(istab0).flat)):

        aprm = azrur[st, :nzbrur[st]]
        bprm = bzrur[st, :nzbrur[st]]
        xzbr = xzbrur[st, :nzbrur[st]]
        szbr = szbrur[st, :nzbrur[st]]

        # If initial sigma, sigz0, solve for pseudo distance, xpskm
        if sigz0 > zerosig:
            jseg = (szbr > sigz0).argmax()
            xpskm = ( sigz0/aprm[jseg] )**(1./bprm[jseg])
        else:
            xpskm = 0
        # Combine the actual and pseudo distances and compute sigz.
        xkm = xpskm + dx
        jseg = (xzbr>xkm).argmax()

        #jseg = np.empty(xkm.shape, dtype=int)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='invalid value encountered in')
            sigz.flat[i]  = sigz0 if dx < 0 else aprm[jseg] * xkm**bprm[jseg]

    return sigz





def pasquil_gifford_y(istab, xkm):
    # sigty.f
    #                th = 0.017453293*(c(istab)-d(istab)*alog(xkm1))
    #               sigy1 = avefac*yrfac*465.11628*xkm1*tan(th)
    # sigty.f
    # c --- Include the ISC parameters for PGT (i.e.,rural) for SIGY.
    #      data c/24.1667,18.333,12.5,8.333,6.25,4.1667/
    #            data d/2.5334,1.8096,1.0857,0.72382,0.54287,0.36191/


    istab0 = istab - 1
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='invalid value encountered in')
        th = (np.pi / 180) * (c[istab0]-d[istab0]*np.log(xkm))
    sigy = 465.11628 * xkm * np.tan(th)
    return sigy


def gaussian(q, ubar, sigy, sigz, h, y, z):
    '''
    q: emission rate (g/sec)
    ubar: wspd (m/sec)
    sigy: (m??)
    sygz: (m??)
    h: height of source (m)
    y: receptor horizontal distance (m)
    z: height of receptor (m)
    '''

    # q, ubar, sig, sigz <= time dependent
    # h <= constant
    # y, z <= space dependent

    c = q / ( 2 * np.pi * ubar * sigy * sigz)
    c *= np.exp( - y*y / ( 2 * sigy*sigy ) )
    c *= ( 
            np.exp( - (z-h)*(z-h) / ( 2 * sigz*sigz ) )
            + 
            np.exp( - (z+h)*(z+h) / ( 2 * sigz*sigz ) )
            )
    return c
