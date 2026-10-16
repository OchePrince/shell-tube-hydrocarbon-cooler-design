"""Hand calculation for the final Design B geometry, using Aspen's fluid properties.

Shell side: Kern method. Tube side: Dittus-Boelter with Coulson-style friction factor.
The results are compared with Aspen EDR in the README (section "Hand calculation versus Aspen EDR").
Run:  python3 hand_vs_edr.py
"""
import math

# Design B geometry (SI units)
Nt = 317                 # number of tubes
do = 0.015875            # tube OD, m (0.625 in)
di = do - 2 * 0.001651   # tube ID, m (0.065 in wall)
Pt = 0.0206375           # triangular pitch, m (0.8125 in)
Ds = 0.457               # shell ID, m (18 in)
B = 0.1651               # baffle spacing, m (6.5 in)
L = 2.438                # tube length, m
Np = 2                   # tube passes

# Hydrocarbon (shell side), Aspen pseudo-component at 150 C
mh, rho, cp, mu, k = 12.0, 753.0, 2181.0, 0.00059, 0.099
# Cooling water (tube side), mean temperature
mc, rc, mu_c, kc, cpc = 28.71, 993.0, 0.00068, 0.628, 4180.0

Q = 1570.3e3             # heat duty, W (Aspen)

# Shell side (Kern)
de = 1.10 / do * (Pt**2 - 0.917 * do**2)       # equivalent diameter, m
As = (Pt - do) * Ds * B / Pt                   # crossflow area, m2
Gs = mh / As
us = Gs / rho
Re = Gs * de / mu
Pr = cp * mu / k
Nu = 0.36 * Re**0.55 * Pr ** (1 / 3)
ho = Nu * k / de

# Tube side
at = Nt / Np * math.pi * di**2 / 4
v = mc / rc / at
Ret = rc * v * di / mu_c
Prt = cpc * mu_c / kc
Nut = 0.023 * Ret**0.8 * Prt**0.4
hi = Nut * kc / di
f = 4 * (0.0014 + 0.125 * Ret**-0.32)
dp_tube = (f * Np * L / di + Np * 2.5) * rc * v**2 / 2


def overall_u(ho_, hi_, rfo=0.0002, rfi=0.00018, kw=50.0):
    """Overall U referred to the tube outside area."""
    return 1 / (1 / ho_ + rfo + do * math.log(do / di) / (2 * kw)
                + do / di * (rfi + 1 / hi_))


Ud = overall_u(ho, hi)
Uc = overall_u(ho, hi, 0.0, 0.0)
lmtd_eff = 80.18                                # effective MTD from EDR, C
A_req = Q / (Ud * lmtd_eff)

print(f"Shell side: velocity {us:.2f} m/s, Re {Re:.0f}, Pr {Pr:.1f}, ho {ho:.0f} W/m2K")
print(f"Tube side:  velocity {v:.2f} m/s, Re {Ret:.0f}, hi {hi:.0f} W/m2K "
      f"({hi*di/do:.0f} on outside area)")
print(f"Tube dP: {dp_tube/1e5:.3f} bar")
print(f"U dirty {Ud:.0f}, U clean {Uc:.0f} W/m2K")
print(f"Required area {A_req:.1f} m2")
