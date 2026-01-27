import CoolProp.CoolProp as cp
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar
import math
import json

with open("evapProps.json") as propsFile:
    props = json.load(propsFile)

# --- global variables (geomerical and thermal properties) ---
nRows = props["No. of Rows"]  # number of rows
nHoles = props["No. of Holes"]  # number of holes
transversePitch = props["Transverse Pitch"] * 1e-3  # transverse tube pitch in m
longitudinalPitch = props["Longitudinal Pitch"] * 1e-3  # longitudinal tube pitch in m
nCircuits = props["No. of Circuits"]  # number of circuits in the evaporator
finThk = props["Fin Thickness"] * 1e-3  # fin thickness in m
tubeThk = props["Tube Thickness"] * 1e-3  # tube thickness
tubeOD = props["Tube OD"] * 1e-3  # tube outer diameter after expansion
finPitch = props["Fin  Pitch"] * 1e-3  # fin pitch
finnedLength = props["Finned Length"] * 1e-3  # finned length

refrigerant = props["Refrigerant"]  # refrigerant
vDot_air = props["Volumetric Air flow"] / 3600.0  # volumetric air flow in m3/sec
atmPressure = props["Atmospheric Pressure"] / 14.5 * 1e5  # atmospheric pressure in Pa
opnPressure = (props["Operating Pressure"] + 14.5) / 14.5 * 1e5  # evaporator pressure
inletAirRH = props["Inlet Air RH"]  # ambient relative humidity (as a fraction)
inletAirTemp = props["Inlet Air Temp"] + 273.15  # ambient temperature in Kelvin

finDepth = nRows * longitudinalPitch  # depth in the air flow direction in m
finHeight = nHoles * transversePitch  # fin height

# calculate the properties of the air
rho_air_in = 1 / cp.HAPropsSI(
    "Vha", "T", inletAirTemp, "P", atmPressure, "RH", inletAirRH
)
mu_air_in = cp.HAPropsSI("M", "T", inletAirTemp, "P", atmPressure, "RH", inletAirRH)

finCollarOD = tubeOD + 2 * finThk  # outer tube diameter including the collars
tubeID = tubeOD - 2 * tubeThk  # internal tube diameter

areaMinCrossSec = (
    finnedLength * (finHeight - nHoles * finCollarOD) * (1 - finThk / finPitch)
)
areaOuter = 2 * (
    finHeight * finDepth - (math.pi * finCollarOD**2 * nHoles * nRows) / 4
) * (finnedLength / finPitch) + (nRows * nHoles * finnedLength) / finPitch * (
    math.pi * (finCollarOD)
) * (
    finPitch - finThk
)  # total outside area = (fin rectangle area - fin holes area)*(number of fins) + area of pipe outer(with fin) * number of pipes
areaInternal = nRows * nHoles * finnedLength * (math.pi * tubeID)
mDot_air = vDot_air * rho_air_in
Gmax_air = mDot_air / areaMinCrossSec
Dh_air = (4 * areaMinCrossSec * finDepth) / areaOuter

Re_air_Dc = (Gmax_air * finCollarOD) / mu_air_in

# --- air-side ---
Pr_air = cp.PropsSI("Prandtl", "T", inletAirTemp, "P", atmPressure, "Air")
Cp_air = cp.PropsSI("C", "T", inletAirTemp, "P", atmPressure, "Air")

# Using Wang et al. (2000) for plain fins
# Colburn j-factor
P3 = (
    -0.361
    - (0.042 * nRows) / math.log(Re_air_Dc)
    + 0.158 * math.log(nRows * (finPitch / finCollarOD) ** (0.41))
)
P4 = -1.224 - (0.0176 * (longitudinalPitch / Dh_air) ** (1.42)) / math.log(Re_air_Dc)
P5 = -0.083 + (0.058 * nRows) / math.log(Re_air_Dc)
P6 = -5.735 + 1.21 * math.log(Re_air_Dc / nRows)
j = (
    0.086
    * Re_air_Dc**P3
    * nRows**P4
    * (finPitch / finCollarOD) ** P5
    * (finPitch / Dh_air) ** P6
    * (finPitch / transversePitch) ** (-0.93)
)

# Friction factor
F1 = (
    -0.764
    + 0.739 * (transversePitch / longitudinalPitch)
    + 0.177 * (finPitch / finCollarOD)
    - 0.00758 / nRows
)
F2 = -15.689 + 64.021 / (math.log(Re_air_Dc))
F3 = 1.696 - 15.695 / math.log(Re_air_Dc)
f_air = (
    0.0267
    * Re_air_Dc**F1
    * (transversePitch / longitudinalPitch) ** F2
    * (finPitch / finCollarOD) ** F3
)

# Calculating heat
h_air = (j * Gmax_air * Cp_air) / (Pr_air ** (2 / 3))
dP_air = f_air * (areaOuter / areaMinCrossSec) * (Gmax_air**2) / (2 * rho_air_in)

# -- refrigerant side ---
nIter = props["Number of elements"]
dl = ((finnedLength * nRows * nHoles) / nCircuits) / nIter
g = 9.81  # accn due to gravity
F_fl = 3.3 # as per Kandlikar (2006)
x_in = props["Inlet Quality"]
mDot_ref = props["Refrigerant Flow Rate"]  # refrigerant mass flow rate [kg/s]
G_ref = mDot_ref / (nCircuits * 0.25 * math.pi * tubeID**2)
P = np.zeros(nIter)  # [evapPressure]
x = np.zeros(nIter)  # [x_in]
H = np.zeros(nIter)  # [cp.PropsSI('H', 'P', evapPressure, 'Q', x_in, refrigerant)]
T = np.zeros(nIter)  # [cp.PropsSI('T', 'P', evapPressure, 'Q', x_in, refrigerant)]
Qtot = np.zeros(nIter)
length = np.zeros(nIter)
localAirTemp = np.full(nIter, inletAirTemp)
mDotPerSegment = mDot_air / (nIter / nRows)

for o in range(5):
    P[0] = opnPressure
    x[0] = x_in
    H[0] = cp.PropsSI("H", "P", opnPressure, "Q", x_in, refrigerant)
    T[0] = cp.PropsSI("T", "P", opnPressure, "Q", x_in, refrigerant)
    length[0] = 0
    Qtot = 0
    QperSegment = np.zeros(nIter)  # We need to store Q to calculate air drop

    for u in range(nIter - 1):
        relaxation = 1
        mu_ref_satV = cp.PropsSI("V", "P", P[u], "Q", 1, refrigerant)
        rho_ref_satV = cp.PropsSI("D", "P", P[u], "Q", 1, refrigerant)
        H_ref_satL = cp.PropsSI("H", "T", T[u], "Q", 0, refrigerant)
        H_ref_satV = cp.PropsSI("H", "T", T[u], "Q", 1, refrigerant)
        latentHeatRef = (
            H_ref_satV - H_ref_satL
        )  # latent heat of vapourisatoin of refrigeration
        Re_satV = (G_ref * tubeID) / mu_ref_satV

        if 0 < x[u] < 1:

            k_satL = cp.PropsSI("L", "T", T[u], "Q", 0, refrigerant)
            mu_satL = cp.PropsSI("V", "P", P[u], "Q", 0, refrigerant)
            Re_satL = (G_ref * tubeID) / mu_satL
            Pr_satL = cp.PropsSI("Prandtl", "T", T[u], "Q", 0, refrigerant)
            rho_ref_satL = cp.PropsSI("D", "P", P[u], "Q", 0, refrigerant)
            rho_ref_homo = 1 / (x[u]/rho_ref_satV + (1-x[u])/rho_ref_satL)
            Fr_lo = (G_ref**2) / (rho_ref_satL**2 * g * tubeID)
            Fr_h  = (G_ref**2) / (rho_ref_homo**2 * g * tubeID)
            h_satL = (
                0.023 * (Re_satL**0.8) * (Pr_satL**0.4) * (k_satL / tubeID)
            )  # Dittus-Boelter
            Co = ((1 - x[u]) / x[u]) ** 0.8 * (
                rho_ref_satV / rho_ref_satL
            ) ** 0.5  # convection number

            def err(T_wall_guess):
                Qguess = (
                    h_air
                    * (areaOuter / (nCircuits * nIter))
                    * (localAirTemp[u] - T_wall_guess)
                )
                qGuess = Qguess / (math.pi * tubeID * dl)
                Bo = qGuess / (G_ref * latentHeatRef)  # boiling number
                h_n = h_satL * (
                    0.6683 * Co ** (-0.2) * (25 * Fr_lo) ** 0.3 + 1058 * Bo**0.7 * F_fl
                )
                h_c = h_satL * (
                    1.136 * Co ** (-0.9) * (25 * Fr_lo) ** 0.3 + 667.2 * Bo**0.7 * F_fl
                )
                href = max(h_c, h_n)
                Qref = href * (math.pi * tubeID * dl) * (T_wall_guess - T[u])
                return Qguess - Qref

            sol = root_scalar(
                err, bracket=[T[u], localAirTemp[u]], method="brentq"
            )  # finding the root to the residual equation using Brent's method
            Twall = sol.root

            Qair = h_air * (areaOuter / (nCircuits * nIter)) * (localAirTemp[u] - Twall)
            heatFlux = Qair / (math.pi * tubeID * dl)
            Bo = heatFlux / (G_ref * latentHeatRef)  # boiling number
            h_nucleate = h_satL * (
                0.6683 * Co ** (-0.2) * (25 * Fr_lo) ** 0.3 + 1058 * Bo**0.7 * F_fl
            )
            h_convective = h_satL * (
                1.136 * Co ** (-0.9) * (25 * Fr_lo) ** 0.3 + 667.2 * Bo**0.7 * F_fl
            )
            h_ref = max(h_convective, h_nucleate)
            Qref = h_ref * (math.pi * tubeID * dl) * (Twall - T[u])

            # test = Qair - Qref
            # print(test)

            H[u + 1] = H[u] + (Qref / mDot_ref)

            f_satL = 0.079 / ((Re_satL) ** 0.25)

            def friedelMultiplier(x):
                f_satV = 0.079 / ((Re_satV) ** 0.25)
                surfTension = cp.PropsSI(
                    "surface_tension", "T", T[u], "Q", 0.5, refrigerant
                )
                WeberNumber = (G_ref**2 * tubeID) / (surfTension * rho_ref_homo)
                E = (1 - x) ** 2 + x**2 * (rho_ref_satL * f_satV) / (
                    rho_ref_satV * f_satL
                )
                Fterm = x**0.78 * (1 - x) ** 0.224
                Hterm = (
                    (rho_ref_satL / rho_ref_satV) ** 0.91
                    * (mu_ref_satV / mu_satL) ** 0.19
                    * (1 - (mu_ref_satV / mu_satL)) ** 0.7
                )
                return E + (3.24 * Fterm * Hterm) / (Fr_h**0.045 * WeberNumber**0.035)

            psi2 = friedelMultiplier(x[u])
            dPLiquidOnly = f_satL * dl / tubeID * G_ref**2 / (2 * rho_ref_satL)
            dP_fric = dPLiquidOnly * psi2

            def voidFraction(x):
                return 1 / (1 + ((1 - x) / x) * (rho_ref_satV / rho_ref_satL))

            x[u + 1] = (H[u + 1] - H_ref_satL) / (H_ref_satV - H_ref_satL)

            def momentumFlux(x, voidFraction):
                term_liq = ((1 - x) ** 2) / (rho_ref_satL * (1 - voidFraction))
                term_vap = (x**2) / (rho_ref_satV * voidFraction)
                return term_liq + term_vap

            Mom_in = momentumFlux(x[u], voidFraction(x[u]))
            Mom_out = momentumFlux(x[u + 1], voidFraction(x[u + 1]))

            dP_mom = G_ref**2 * (Mom_out - Mom_in)

            dP = dP_fric + dP_mom

        else:
            mu_v_sh = cp.PropsSI("V", "P", P[u], "H", H[u], refrigerant)
            rho_v_sh = cp.PropsSI("D", "P", P[u], "H", H[u], refrigerant)
            Re_v_sh = (G_ref * tubeID) / mu_v_sh
            Pr_v_sh = cp.PropsSI("Prandtl", "H", H[u], "P", P[u], refrigerant)
            k_v_sh = cp.PropsSI("L", "H", H[u], "P", P[u], refrigerant)
            if Re_v_sh < 3000:
                Nu_v_sh = 4.36
                f_ref_v = 64 / Re_v_sh
            else:
                f_ref_v = (0.790 * math.log(Re_v_sh) - 1.64) ** (-2)
                Nu_v_sh = ((f_ref_v / 8) * (Re_v_sh - 1000) * Pr_v_sh) / (
                    1 + 12.7 * (f_ref_v / 8) ** (1 / 2) * (Pr_v_sh ** (2 / 3) - 1)
                )

            h_ref = Nu_v_sh * k_v_sh / tubeID

            if T[u] >= localAirTemp[u] - 0.001:  # 0.001 tolerance for float precision
                Twall = localAirTemp[u]
                Qref = 0
            else:
                # Only run solver if there is a temperature difference
                def err(T_wall_guess):
                    Qguess = (
                        h_air
                        * (areaOuter / (nCircuits * nIter))
                        * (localAirTemp[u] - T_wall_guess)
                    )
                    Qref = h_ref * (math.pi * tubeID * dl) * (T_wall_guess - T[u])
                    return Qguess - Qref

                try:
                    # Use a slightly padded upper bound to avoid exact edge cases
                    sol = root_scalar(
                        err, bracket=[T[u], localAirTemp[u]], method="brentq"
                    )
                    Twall = sol.root
                    Qref = h_ref * (math.pi * tubeID * dl) * (Twall - T[u])
                except ValueError:
                    # Fallback if solver fails despite check (rare numerical noise)
                    Twall = T[u]
                    Qref = 0

            H[u + 1] = H[u] + (Qref / mDot_ref)

            dP = f_ref_v * (dl / tubeID) * (G_ref**2 / (2 * rho_v_sh))

            x[u + 1] = 1

        P[u + 1] = P[u] - dP  # *relaxation

        T[u + 1] = cp.PropsSI("T", "P", P[u + 1], "H", H[u + 1], refrigerant)

        QperSegment[u] += Qref
        Qtot += Qref
        length[u + 1] = length[u] + dl

        # if not 0 < x[u+1] < 1:
        #     print('Single Phase Start')
        #     break
    print(f"Pass {o+1}: Total Capacity = {Qtot:.2f} W")

    step = nIter // nRows
    localAirTemp[(nRows - 1) * step :] = inletAirTemp

    for p in range(nRows - 1, 0, -1):
        upwind_start = p * step
        upwind_end = (p + 1) * step

        downwind_start = (p - 1) * step
        downwind_end = p * step

        for k in range(step):

            idx_upwind = upwind_start + k
            idx_downwind = downwind_start + k

            Cp_eff = 2000.0 if inletAirRH > 0.4 else 1006.0

            dT_air = QperSegment[idx_upwind] / (mDotPerSegment * Cp_eff)

            localAirTemp[idx_downwind] = localAirTemp[idx_upwind] - dT_air


print(f"Final Total Capacity = {Qtot:.2f} W")


# Generate saturation lines for plotting
Psat = np.linspace(10e5, 5e5, 500)
h_f_sat = [cp.PropsSI("H", "P", P, "Q", 0, refrigerant) for P in Psat]
h_g_sat = [cp.PropsSI("H", "P", P, "Q", 1, refrigerant) for P in Psat]

h_f_sat = np.array(h_f_sat)
h_g_sat = np.array(h_g_sat)

# Plotting the P-h diagram
plt.figure(figsize=(10, 8))
# Plot saturation lines
plt.plot(h_f_sat * 1e-3, Psat, "b-", label="Saturated Liquid")
plt.plot(h_g_sat * 1e-3, Psat, "r-", label="Saturated Vapor")
plt.plot(H * 1e-3, P)
plt.xlabel("Enthalpy (KJ/kg)")
plt.ylabel("Pressure (Pa)")
plt.title("P-h Diagram of the Evaporation process")
plt.minorticks_on() # 1. Enable minor ticks
plt.grid(which='major', linestyle='-', linewidth='0.5', color='black') # 2. Style major grid
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')  # 3. Style minor grid
plt.savefig("Ph.png", dpi = 300)

plt.figure(figsize=(10, 8))
plt.plot(length, x)
plt.xlabel("Length along the tube (m)")
plt.ylabel("Vapour Quality")
plt.grid(True)
plt.title("Vapour Quality along the length of the tube")
plt.minorticks_on() # 1. Enable minor ticks
plt.grid(which='major', linestyle='-', linewidth='0.5', color='black') # 2. Style major grid
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')  # 3. Style minor grid
plt.savefig("x.png", dpi = 300)

plt.figure(figsize=(10, 8))
plt.plot(length, T - 273.15)
plt.xlabel("Length along the tube (m)")
plt.ylabel("Temperature (degree C)")
plt.grid(True)
plt.title("Temperature profile along the length of the tube")
plt.minorticks_on() # 1. Enable minor ticks
plt.grid(which='major', linestyle='-', linewidth='0.5', color='black') # 2. Style major grid
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')  # 3. Style minor grid
plt.savefig("tempProfile.png", dpi =300)
