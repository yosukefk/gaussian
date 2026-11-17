import numpy as np

def pgtstb(u10, v10, icc=0, iceil=np.inf, sinalp=1, iwat=0):
    '''
    PURPOSE:  Calculate PGT stability class (1-6)

    u10    - real    - Surface U component wind (m/s) at height of 10 m 
    u10    - real    - Surface V component wind (m/s) at height of 10 m 
    ICC    - integer - Cloud cover (tenths) 
    ICEIL  - integer - Ceiling height (hundreds of ft)
    SINALP - real    - Sine of the solar elevation angle
    iwat   - integer - Water flag (1 over water, 0 over land)
    '''
    scrit = np.array([0.258819,0.5735764,0.8660254])
    jstab = np.array([
                 6,6,6,6,6,6,5,5,5,5,4,4,
                 6,6,6,5,5,5,4,4,4,4,4,4,
                 4,4,4,4,4,4,4,4,4,4,4,4,
                 3,3,3,4,4,4,4,4,4,4,4,4,
                 2,2,2,3,3,3,3,3,3,4,4,4,
                 1,2,2,2,2,2,2,3,3,3,3,4,
                 1,1,1,1,1,2,2,2,2,3,3,3,
                 ]).reshape((7,12))



    if iwat != 0: return np.nan

    if sinalp < 0:
        # nightime stability

        # for icc = 0-4, ic = -2
        ic = -2

        if icc > 4: 
            # for icc = 4-9 or 10 (with iceil ge 7,000 ft), icc = -1
            ic = -1 

        if icc == 10 and iceil < 70: 
            # for icc = 10 and iceil lt 7,000 ft, icc = 0
            ic = 0
    else: 
        # daytime stability

        # calculate radiation index, iri 
        # scrit(1)=sin(15. deg)
        # scrit(2)=sin(35. deg)
        # scrit(3)=sin(60. deg)
        iri = 3

        for n in range(3):
            if sinalp < scrit[n]:
                iri = n

# ---    determine daytime insolation class, ic
        ic = iri + 1

        if icc > 5:
            if icc == 10:
# ---          icc = 10  --  ceiling ht. is in hundreds of ft.
                if iceil < 70:
                    # icc=10; iceil lt 7,000 ft.  
                    ic = 0
                elif iceil < 160:
                    # icc=10; 7,000 ft le iceil lt 16,000 ft
                    ic = iri - 2
                    ic = max(ic, 1)
                else:
                    # icc=10; iceil ge 16,000 ft 
                    ic = iri - 1 
                    ic = max(ic, 1)
            else:
                # 5 lt icc lt 10
                if iceil < 70: 
                    # 5 lt icc lt 10; iceil lt 7,000 ft 
                    ic = iri - 2 
                    ic = max(ic, 1)
                elif icel < 160:
                    # 5 lt icc lt 10; 7,000 ft le iceil lt 16,000 ft 
                    ic = iri - 1 
                    ic = max(ic, 1)
                else:
                    # 5 lt icc lt 10; iceil ge 16,000 ft 
                    ic = iri
    # --- convert insolation class (-2 to +4) to array index (1 to 7)
    idx = ic + 3

    # convert wind speed from m/s to nearest knot(1.9438 knots=1 m/s) 
    iwsx = int(1.9438 * np.sqrt(u10*u10+v10*v10) + 0.5 )
    iwsx = max(iwsx,1) 
    iwsx = min(iwsx,12) 
    
    # lookup :  -1 for python indexing
    ipgt = jstab[idx-1, iwsx-1]
    return ipgt





