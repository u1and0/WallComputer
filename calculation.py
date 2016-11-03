"""
## calculation ver1.0

__USAGE__
Just build

__INTRODUCTION__

'Analysis of Radome-enclosed Antennas'
Author : D. J. Kozakoff

P85~
'Program Wall Computer Software Listing'

This is copy of Fortran program.


__ACTION__
--action--

__UPDATE1.0__
First commit

__TODO__
None
"""


import numpy as np
# __GLOBAL ARGUMENTS__________________________
Amag, Bmag, Cmag, Aph, Bph, Cph = np.zeros(10), np.zeros(
    10), np.zeros(10), np.zeros(10), np.zeros(10), np.zeros(10)
zmag, zph, rez, imz, emag, eph, er, ltan = np.zeros(10), np.zeros(
    10), np.zeros(10), np.zeros(10), np.zeros(10), np.zeros(10), np.zeros(10)

Twdb = np.zeros((179, 90))  # transmission coefficient(dB)
IPD = np.zeros((179, 90))
N = 10  # Number of layer


def rect_to_polar(x, y):
    mag = np.sqrt(x**2 + y**2)
    corr = 0 if x >= 0 else np.pi
    ang = np.arctan2(y, x) + corr if x != 0 and y != 0 else 0
    return (mag, ang)


def complex_sqr_root(ang):
    rm = np.sqrt(ang)
    ra = ang / 2
    return (rm, ra)


def polar_to_rect(ang, mag):
    x = mag * np.cos(ang)
    y = mag * np.sin(ang)
    return (x, y)

def cmult():
    for i in range(1,3,2)


def wall(theta):
    CTH = np.cos(theta)
    STH = np.sin(theta)**2
    zmag[0], zph[0], rez[0], imz[0] = 1, 0, 1, 0
    # __Calculate Complex Permitivity For Each Layer__________________________
    for I in range(1, N + 1):
        x = er[I]
        y = er[I] * -ltan[I]
        MAG, ANG = rect_to_polar(x, y)
        emag[I] = MAG
        eph[I] = ANG
        RMAG, RANG = complex_sqr_root(ANG)
        magterm[I], angterm[I] = RMAG, RANG
        phimag[I] = Ko * magterm[I]
        phiph[I] = angterm[I]
        MAG, ANG = phimag[I], phiph[I]
        x, y = polar_to_rect()
        rephi[I] = x
        imphi[I] = y
        zmag[I] = CTH / magterm[I]
        x, y = polar_to_rect(ANG, MAG)
        rez[I], imz[I] = x, y
        if pol == 1:
            zmag[I] = 1 / (emag[I], zmag[I])
            zph[I] = -(eph[I] + zph[I])
            MAG, ANG = zmag[I], zphi[I]
            x, y = polar_to_rect(ANG, MAG)
            rez[I], imz[I] = x, y
    zmag[N + 1], zph[N + 1], rez[N + 1], imz[N + 1] = 1, 0, 1, 0
    for I in range(1, N + 2):
        """最適化前
        x = rez[I] - rez[I - 1]
        y = imz[I] - imz[I - 1]
        MAG, ANG = rect_to_polar(x, y)
        NUMMAG, NUMANG = MAG, ANG
        x = rez[I] + rez[I - 1]
        y = imz[I] + imz[I - 1]
        MAG, ANG = rect_to_polar(x, y)
        DENMAG, DENANG = MAG, ANG
        RMAG[I] = NUMMAG / DENANG
        Rph[I] = NUMMAG - DENANG
        MAG, ANG = RMAG[I], Rph[I]
        x, y = polar_to_rect(ANG, MAG)
        reR[I], imR[I] = x, y
        reT[I], imT[I] = 1 + reR[I], imR[I]
        x, y = reT[I], imT[I]
        MAG, ANG = rect_to_polar(x, y)
        Tmag[I], Tph[I] = MAG, ANG"""
        NUMMAG, NUMANG = rect_to_polar(rez[I] - rez[I - 1], imz[I] - imz[I - 1])
        DENMAG, DENANG = rect_to_polar(rez[I] + rez[I - 1], imz[I] + imz[I - 1])
        RMAG[I] = NUMMAG / DENANG
        Rph[I] = NUMMAG - DENANG
        reR[I], imR[I] = polar_to_rect(RMAG[I], Rph[I])
        reT[I], imT[I] = 1 + reR[I], imR[I]
        Tmag[I], Tph[I] = rect_to_polar(reT[I], imT[I])
    Amag[1] = np.exp(-imphi[1] * thk[1])
    Amag[4] = 1 / Amag[1]
    Amag[2] = RMAG[1] * Amag[4]
    Amag[3] = RMAG[1] * Amag[1]
    Aph[1] = ephi[1] * thk[1]
    Aph[2] = Rph[1] - Aph[1]
    Aph[3] = Rph[1] + Aph[1]
    Aph[4] = -Aph[1]
    for K in range(2, N + 1):
        Bmag[1] = np.exp(-imphi[K] * thk[K])
        Bmag[4] = 1 / Bmag[1]
        Bmag[2] = RMAG[K] * Bmag[4]
        Bmag[3] = RMAG[K] * Bmag[1]
        Bph[1] = ephi[K] * thk[K]
        Bph[2] = Rph[K] - Bph[1]
        Bph[3] = Rph[K] + Bph[1]
        Bph[4] = -Bph[1]
        cmult()
# __comment...__________________________


def main():
    # __CALCULATION START__________________________
    array = np.array([])
    """fl, fh, finc
    あとでinput()させる。
    とりあえず値を入れておく。
    単位はGHz
    """
    fl, fh, finc = 1, 10, 1
    for pol in range(1):
        array[0] = 0
        for f in np.arange(fl, fh, finc):
            array += 1

            LAM = 11.811 / f
            Ko = .531976 * f

        for ANGLE in np.arange(0., 90., 15.):
            theta = np.radians(ANGLE)
            wall(theta)


if __name__ == '__main__':
    main()
