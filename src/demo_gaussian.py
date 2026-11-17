from gaussian import gaussian as g
import numpy as np


# %%%%%%%%%%%%%%
# sample data
# %%%%%%%%%%%%%%

ws = 3 # m/sec
#wd = 90 # from east

#wd = 0 # from north
#wd = 45 # from NE
#wd = 315 # from NW
wd =  225 # from SW

# coordinates, in meters
src = np.array((0, 0))
rcp = np.array((50, 30))

# stability class
istab = 3

# emission rate
emis = 10 # kg/hr

# height of souce
h = 5.5

# height of receoptor
z = 2

print(f'''
Sample Data:
ws = {ws} m/sec
wd = {wd} deg
istab = {istab}

src coords = {src}
rcp coords = {rcp}
h = {h} m (source)
z = {z} m (receptor)

emis = {emis} kg/hr
'''
)



# get angle
dst = np.sqrt(((src - rcp ) ** 2).sum())

def wd2rad(wd):
    # degree in wd
    # returns unit vector representing the wind
    wrad = (270 - wd) / 180 * np.pi
    return wrad

def xy2rad(x, y):
    # xy displacement vector
    # returns angle
    drad = np.arctan2(x, y)
    return drad

wrad = wd2rad(wd)
drad = xy2rad( (rcp-src)[0], (rcp-src)[1])

rrad = drad - wrad

print(f'''
angles in radians:  
    wind vector        = {wrad:.3f} rad
    directional vector = {drad:.3f} rad
    rotational vector  = {rrad:.3f} rad
angles in degrees:
    wind vector        = {wrad*180/np.pi:.1f} deg
    directional vector = {drad*180/np.pi:.1f} deg
    rotational vector  = {rrad*180/np.pi:.1f} deg
    '''
    )

# how should i call longitudinal/lateral distnce....?
u = dst * np.cos(rrad)
v = dst * np.sin(rrad)
print(f'''
distance:
    longitudinal = {u:.1f} m
    latetal      = {v:.1f} m
    ''')



sigy = g.pasquil_gifford_y(istab, u/1000.)
sigz = g.pasquil_gifford_z(istab, u/1000.)

#print(f'sigy={sigy}, sigz={sigz}')
print(f'''
dispersion coeffs:
    sigy = {sigy:.2f}
    sigz = {sigz:.2f}
    ''')

q = emis * 1000 / 3600.

conc = g.gaussian(q, ws,  sigy, sigz, h, v, z)
# convert unit of array from g/m3 tp ppb
# mwt g/mol
# molar volume m3/mol
conc = conc / 16.043 * 0.024465403697038 * 1e9
#print(f'conc={conc}')
print(f'''
predicted conc:
    conc = {conc:.1f} ppb
    ''')
