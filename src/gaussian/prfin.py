import numpy as np
def prfin(vexit,diam,fm,fb,wsstk,istab,sqrts,):
    ''' 
    PURPOSE:  Calculate distance to final rise and height of final plume rise
    INPUTS:
          VEXIT - real         - Stack gas exit velocity (m/s)
           DIAM - real         - Stack diameter (m)
             FM - real         - Momentum flux (m**4/s**2)
             FB - real         - Buoyancy flux (m**4/s**3)
          WSSTK - real         - Stack height wind speed (m/s)
          ISTAB - integer      - Stability class (1=A,...6=F)
          SQRTS - real         - Square root of stability parameter,s
                            (sqrt(s)=(g*(dtheta/dz)/tair)**0.5)
      OUTPUTS:
            XFINM - real         - Distance (m) to final momentum
                                   plume rise
            XFINB - real         - Distance (m) to final buoyant
                                   plume rise
             XFIN - real         - Larger of XFINM and XFINB
            ZFINM - real         - Final momentum plume rise (m)
            ZFINB - real         - Final buoyant plume rise (m)
             ZFIN - real         - Larger of ZFINM and ZFINB
    '''

# --- Set minimum allowed wind speed at stack-top (wsstk0).  If wind
# --- speed is less than wsstk0, set distance to final rise equal to
# --- zero, and use the "calm rise" formula for stable conditions, or
# --- the neutral/unstable rise formulas with wsstk=wsstk0.
    wsstk0 = 1
# ------------------------------------------
# --- Compute distance, height of final rise
# ------------------------------------------
    if wsstk >= wsstk0: 
# ---    "non-CALM" conditions ------------------------

# ---    Buoyant plume rise
        if istab <= 4: 
# ---       Neutral/unstable conditions
            if fb >= 55.: 
                xfinb = 119. * (fb ** 0.4) 
                zfinb = 38.71 * (fb ** 0.6) / wsstk
            elif fb > 0.0:
                xfinb = 49. * (fb ** 0.625 )
                zfinb = 21.425 * (fb ** 0.75) / wsstk
            else: 
# ---          XFINB with FB le 0.0 follows ISC2 convention
                vexit2 = max(vexit,1.0e-10) 
                xfinb = 4. * diam * ( ( vexit2 + 3.*wsstk) ** 2) / (vexit2*wsstk) 
                fb2 = max(fb,1.0e-10) 
                zfinb = 21.425 * (fb ** 0.75) / wsstk
        else:
#
# ---       Stable conditions
            xfinb = 2.0715 * wsstk / sqrts
            s = sqrts * sqrts
            fb2 = max(fb,1.0e-10)
#
# ---       Bent-over plume equation
            zfinb = 2.6 * ( fb2 / (wsstk*s) ) ** 0.333333
#
# ---       Vertical plume equation
            zfinv = 4. * (fb2 ** 0.25) / (s ** 0.375)
#
# ---       Final stable rise is lower of bent-over and vertical plume
# ---       values
            zfinb=min(zfinb,zfinv)

#
# ---    Momentum plume rise
        if istab <= 4: 
#
# ---       Neutral/unstable conditions
            vexit2 = max(vexit,1.0e-10)
            xfinm = 4. * diam * ((vexit2+3.*wsstk)**2)/(vexit2*wsstk)
            zfinm = 3. * diam * vexit/wsstk
        else: 
#
# ---       Stable conditions
# ---       (Constant 1.57080 = pi/2.)
            xfinm=1.57080*wsstk/sqrts
            fm2=max(fm,1.0e-10)
            zfinm=1.5*(fm2/(wsstk*sqrts))**0.333333
#
# ---       Maximum stable momentum rise is restricted to be less than
# ---       final neutral/unstable momentum rise
            zfinm=min(zfinm,3.*diam*vexit/wsstk) 
# ---    Set xfinm,zfinm to zero explicitly when FM is zero
        if fm == 0.0:
            xfinm=0.0
            zfinm=0.0

    else: 
# ---    "CALM" conditions ----------------------------
#
# ---    Buoyant plume rise
        if istab <= 4: 
#
# ---       Neutral/unstable conditions
            xfinb=0.0
            if fb >= 55.: 
                zfinb=38.71*(fb**0.6)/wsstk0 
            else: 
                fb2=max(fb,1.0e-10) 
                zfinb=21.425*(fb2**0.75)/wsstk0 
        else:
#
# ---       Stable conditions
            xfinb=0.0
            s=sqrts*sqrts
            fb2=max(fb,1.0e-10)
#
# ---       Vertical plume equation
            zfinb=4.*(fb2**0.25)/(s**0.375) 
#
# ---    Momentum plume rise
        if istab >= 4: 
#
# ---       Neutral/unstable conditions
            vexit2=max(vexit,1.0e-10)
            xfinm=0.0
            zfinm=3.*diam*vexit/wsstk0
        else:
#
# ---       Stable conditions
            xfinm=0.0
            fm2=max(fm,1.0e-10)
            zfinm=1.5*(fm2/(wsstk0*sqrts))**0.333333
#
# ---       Maximum stable momentum rise is restricted to be less than
# ---       final neutral/unstable momentum rise
            zfinm=min(zfinm,3.*diam*vexit/wsstk0)
# ---    Set xfinm,zfinm to zero explicitly when FM is zero 
        if fm == 0.0: 
            xfinm=0.0
            zfinm=0.0
#
# --- Set overall maximum distance to final plume rise and maximum
# --- plume rise 
    xfin=max(xfinb,xfinm) 
    zfin=max(zfinb,zfinm)
    return xfinm, xfinb, xfin, zfinm, zfinb, zfin

def briggs_rise(diam, exitw,tstak,temp, ws, istab): 
    '''
    DIAM - real     - Stack diameter (m)
    EXITW - real     - Exit velocity (m/s)
    TSTAK - real     - Exit temperature (deg. K)
    TEMP  -  2m tempr?  or at stack tip?(deg.K)
    WS    - 10m tempr?  or at stack tip? m/sec
    ISTAB stability class
    '''

    # SQRTS calculations
    # from rlsmet.f
    g = 9.80665
# --- Set minimum allowed potential temperature gradient PTGRAD0,
# --- and the layer above stack-top over which it is calculated DELZ
    ptgrad0 = .005
    # from  calpuff.inp
#        Default potential temperature gradient
#        for stable classes E, F (degK/m)
#        (PTG0(2))                       Default: 0.020, 0.035
    ptg0 = [0.020, 0.035]
    if istab <= 4:
        ptgrad = 0.0
    elif istab ==5:
        ptgrad = max(ptg0[0], prgrad0)
    else:
        ptgrad = max(ptg0[1], prgrad0)
    sqrts=np.sqrt(g*ptgrad/temp)

    # 
    # Set temperature excess at release
    dtemp=tstak-temp
    dtemp=max(dtemp,0.0)

    # Calculate fluxes
    # Buoyancy flux:  2.4516625 = g/4., where g = 9.80665 m/s**2
    fluxb=(2.4516625*diam**2*exitw/tstak)*dtemp
    # Momentum flux 
    fluxm=0.25*(exitw*diam)**2*temp/tstak
    # --- Final plume rise (no downwash considered)
    (xfinm_B,xfinb_B,xf_B,zfinm_B,zfinb_B,rise_B) = prfin(exitw, diam, fluxm, fluxb, ws, istab, sqrts,)

    # briggs rise
    rise = rise_B
    return rise
