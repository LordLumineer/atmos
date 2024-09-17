"""
Find gas properties in earth's atmosphere.

This module provides a function `atmos` which computes the gas properties in
earth's atmosphere. The function takes as input the geopotential altitude in
meters and a temperature offset in degrees Celsius. It can also take as option
the unit system for inputs and outputs and the use of geometric instead of
geopotential altitude.

The function returns a tuple containing the density, speed of sound,
temperature, pressure, kinematic viscosity and height or altitude of the
atmosphere at the given input.

MIT License

Copyright (c) 2024 LordLumineer

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import numpy as np
import pint

ureg = pint.UnitRegistry(cache_folder=":auto:", system='SI')
Q_ = ureg.Quantity


def atmos(
    H_in: Q_ | float = Q_(0, ureg.meter),
    T_offset: Q_ | float = Q_(0, ureg.delta_degC),
    Units: list[str] | str = 'SI',
    GeomFlag: bool = False
) -> list:
    """
    ATMOS Find gas properties in earth's atmosphere.

    Parameters
    ----------
    - H_in : float or pint.Quantity[length], optional
        - Geopotential altitude in meters or feet. 
        - Default is 0 meters.
    - T_offset : float or pint.Quantity[delta_Temp], optional
        - Temperature offset in degrees Celsius or Fahrenheit (or absolute Kelvin or Rankine). 
        - Default is 0.
    - Units : str or list[str], optional
        - Unit system for inputs and outputs. 
        - Default is 'SI'.
    - GeomFlag : bool, optional
        - If True, use geometric altitude instead of geopotential altitude. 
        - Default is False.

    Returns
    -------
    rho : float or pint.Quantity
        Density in kilograms per cubic meter or slugs per cubic foot.
    a : float or pint.Quantity
        Speed of sound in meters per second or feet per second.
    T : float or pint.Quantity
        Temperature in degrees Kelvin or Rankine.
    P : float or pint.Quantity
        Pressure in pascals or pounds per square foot.
    kvisc : float or pint.Quantity
        Kinematic viscosity in square meters per second or square feet per second.
    ZorH : float or pint.Quantity
        Height or altitude in meters or feet.

    Example
    -------
    >>> from atmos import atmos
    >>> rho, a, T, P, kvisc, ZorH = atmos(H_in=1000, T_offset=0, Units='SI', GeomFlag=False)
    >>> print("Density:", rho)
    >>> print("Speed of sound:", a)
    >>> print("Temperature:", T)
    >>> print("Pressure:", P)
    >>> print("Kinematic viscosity:", kvisc)
    >>> print("Height or altitude:", ZorH)

    Notes
    -----
    ATMOS Find gas properties in earth's atmosphere.

    ATMOS by itself gives the atmospheric properties at sea level on a standard day.

    ATMOS(H) returns the properties of the 1976 Standard Atmosphere at geopotential altitude H (meters).

    ATMOS(H,dT) returns properties when the temperature is dT degrees offset from standard conditions.

    ATMOS(H,dT,Units) specifies units for the inputs outputs.

    - The unit system can be specified as a string or a list of strings. 
    - If a single string is provided, it will be used for both the input and output units.  
    - If a list of two strings is provided, the first string will be used for the input units and the second string will be used for the output units.  
    - The supported unit systems are: SI (default), US (a.k.a. Imperial, English), (others pint supported unit systems are not yet supported).  
    - If a pint.Quantity is provided as an input, the units will be automatically detected.  

    For example:
        To use the SI system, set Units to '', or 'SI'.
        To use the US system, set Units to 'US'.
        To use different units for the input and output, pass a list of strings, e.g. ['US', 'SI'].

    The temperature offset (dT) is in the same units as the input temperature, 
    so when converting between Celsius and Fahrenheit, use only the scaling factor (dC/dF = dK/dR = 5/9).

    Units are as follows:
        Input:                        SI (default)      US
            H:     Altitude           m                 ft
            dT:    Temp. offset       K/delta_degC      R/delta_degF
        Output:
            rho:   Density            kg/m^3            slug/ft^3
            a:     Speed of sound     m/s               ft/s
            T:     Temperature        degK              degR
            P:     Pressure           Pa                lbf/ft^2
            nu:    Kinem. viscosity   m^2/s             ft^2/s
            ZorH:  Height or altitude m                 ft

    This atmospheric model is not recommended for use at altitudes above
    86 km geometric height (84852 m/278386 ft geopotential) and returns NaN
    for altitudes above 100 km geopotential.

    References: ESDU 77022;  www.pdas.com/atmos.html
    """

    # Check Units
    if isinstance(Units, str):
        # Same Units for input and output
        Units = 'SI' if Units == '' else Units
        # dir(ureg.sys) not yet supported
        assert Units.upper() in ['SI', 'US',
                                 ''], "Units must be 'SI', 'US', or '' (-> SI)"
        Units = [Units.upper(), Units.upper()]

    elif isinstance(Units, list):
        if len(Units) == 1:
            # Same Units for input and output
            Units = Units[0]
            Units = 'SI' if Units == '' else Units
            # dir(ureg.sys) not yet supported
            assert Units.upper() in [
                'SI', 'US', ''], "Units must be 'SI', 'US', or '' (-> SI)"
            Units = [Units.upper(), Units.upper()]
        elif len(Units) == 2:
            # Different Units for input and output
            assert Units[0].upper() in [
                # dir(ureg.sys) not yet supported
                'SI', 'US', ''], "Units[0] must be 'SI', 'US', or '' (-> SI)"
            # dir(ureg.sys) not yet supported
            assert Units[1].upper() in ['SI', 'US',
                                        ''], "Units[1] must be 'SI', 'US', or '' (-> SI)"
            Units = [Units[0].upper(), Units[1].upper()]
        else:
            raise ValueError(
                'Units must be a list of length 1 or 2 (Input/Output)')

    else:
        raise ValueError(
            'Units must be a list of length 1 or 2 (Input/Output)')

    # Unit Conversion IN
    dimensioned = False
    if isinstance(H_in, Q_):
        dimensioned = True
        if H_in.check('[length]'):
            H_in = H_in.to(ureg.meter).magnitude
        else:
            raise ValueError('H_in must be a pint Quantity of length')
    elif isinstance(H_in, (int, float)):
        if Units[0] == 'US':
            H_in = Q_(H_in, ureg.ft).to(ureg.meter).magnitude
        else:
            H_in = Q_(H_in, ureg.meter).magnitude
    else:
        raise ValueError('H_in must be a pint Quantity of length or a number')

    if isinstance(T_offset, Q_):
        if T_offset.check('[temperature]'):
            dimensioned = True
            T_offset = T_offset.to(ureg.delta_degC).magnitude
        else:
            raise ValueError('T_offset must be a pint Quantity of temperature')
    elif isinstance(T_offset, (int, float)):
        if Units[0] == 'US':
            T_offset = Q_(T_offset, ureg.delta_degF).to(
                ureg.delta_degC).magnitude
        else:
            T_offset = Q_(T_offset, ureg.delta_degC).magnitude

    else:
        raise ValueError(
            'T_offset must be a pint Quantity of temperature or a number')

    # Constants
    R = 287.05287   # N-m/kg-K; value from ESDU 77022
    gamma = 1.4
    g0 = 9.80665    # m/sec^2
    RE = 6356766    # Radius of the Earth, m
    Bs = 1.458e-6   # N-s/m^2 K^1/2
    S = 110.4       # K

    # Quick troposphere-only code
    if (H_in <= 11000) and (T_offset == 0 and Units[0] == 'SI' and not GeomFlag):
        K_0, T_0, H_0, P_0 = -0.0065, 288.15, 0, 101325

        TonTi = 1 + K_0*(H_in - H_0) / T_0
        temp = TonTi * T_0
        PonPi = TonTi**(-g0 / (K_0 * R))
        press = P_0 * PonPi

        rho = press/(temp*R)
        a = np.sqrt(gamma*R*temp)
        kvisc = (Bs * temp**1.5 / (temp + S)) / rho
        ZorH = RE*H_in/(RE-H_in)

        if dimensioned:
            rho = Q_(rho, ureg.kg / ureg.m**3)
            a = Q_(a, ureg.m / ureg.s)
            temp = Q_(temp, ureg.K)
            press = Q_(press, ureg.Pa)
            kvisc = Q_(kvisc, ureg.m**2 / ureg.s)
            ZorH = Q_(ZorH, ureg.m)
            if Units[1] == 'US':
                rho = rho.to(ureg.slug / ureg.ft**3)
                a = a.to(ureg.foot / ureg.second)
                temp = temp.to(ureg.degF)
                press = press.to(ureg.psi)
                kvisc = kvisc.to(ureg.foot**2 / ureg.second)
                ZorH = ZorH.to(ureg.foot)
        return rho, a, temp, press, kvisc, ZorH

    # Lapse rate  Base Temp   Base Geopo Alt    Base Pressure
    # Ki (C/m)    Ti (K)      Hi (m)            P (Pa)
    D = np.array([
        [-0.0065, 288.15, 0, 101325],
        [0, 216.65, 11000, 22632.0400950078],
        [0.001, 216.65, 20000, 5474.87742428105],
        [0.0028, 228.65, 32000, 868.015776620216],
        [0, 270.65, 47000, 110.90577336731],
        [-0.0028, 270.65, 51000, 66.9385281211797],
        [-0.002, 214.65, 71000, 3.9563921603966],
        [0, 186.94590831019, 84852.0458449057, 0.373377173762337]
    ])

    K = D[:, 0]  # K/m
    T = D[:, 1]  # K
    H = D[:, 2]  # m
    P = D[:, 3]  # Pa

    temp = np.zeros_like(H_in)
    press = np.zeros_like(H_in)
    hmax = 100000

    if GeomFlag:
        Hgeop = (RE * H_in) / (RE + H_in)
    else:
        Hgeop = H_in

    n1 = Hgeop <= H[1]
    n2 = Hgeop <= H[2] and Hgeop > H[1]
    n3 = Hgeop <= H[3] and Hgeop > H[2]
    n4 = Hgeop <= H[4] and Hgeop > H[3]
    n5 = Hgeop <= H[5] and Hgeop > H[4]
    n6 = Hgeop <= H[6] and Hgeop > H[5]
    n7 = Hgeop <= H[7] and Hgeop > H[6]
    n8 = Hgeop <= hmax and Hgeop > H[7]
    n9 = Hgeop > hmax

    # Troposphere
    if n1:
        TonTi = 1 + K[0]*(Hgeop - H[0]) / T[0]
        temp = TonTi * T[0]
        PonPi = TonTi**(-g0 / (K[0] * R))
        press = P[0] * PonPi

    # Tropopause
    if n2:
        temp[n2] = T[1]
        PonPi = np.exp(-g0 * (Hgeop - H[1]) / (T[1] * R))
        press = P[1] * PonPi

    # Stratosphere 1
    if n3:
        TonTi = 1 + K[2] * (Hgeop - H[2]) / T[2]
        temp = TonTi * T[2]
        PonPi = TonTi**(-g0 / (K[2] * R))
        press = P[2] * PonPi

    # Stratosphere 2
    if n4:
        TonTi = 1 + K[3] * (Hgeop - H[3]) / T[3]
        temp = TonTi * T[3]
        PonPi = TonTi**(-g0 / (K[3] * R))
        press = P[3] * PonPi

    # Stratopause
    if n5:
        temp = T[4]
        PonPi = np.exp(-g0 * (Hgeop - H[4]) / (T[4] * R))
        press = P[4] * PonPi

    # Mesosphere 1
    if n6:
        TonTi = 1 + K[5] * (Hgeop - H[5]) / T[5]
        temp = TonTi * T[5]
        PonPi = TonTi**(-g0 / (K[5] * R))
        press = P[5] * PonPi

    # Mesosphere 2
    if n7:
        TonTi = 1 + K[6] * (Hgeop - H[6]) / T[6]
        temp = TonTi * T[6]
        PonPi = TonTi**(-g0 / (K[6] * R))
        press = P[6] * PonPi

    # Mesopause
    if n8:
        temp = T[7]
        PonPi = np.exp(-g0 * (Hgeop - H[7]) / (T[7] * R))
        press = P[7] * PonPi

    if n9:
        raise ValueError('altitude too high')

    temp = temp + T_offset
    rho = press / (R * temp)
    a = np.sqrt(gamma * R * temp)
    kvisc = (Bs * temp**1.5 / (temp + S)) / rho
    if GeomFlag:
        ZorH = Hgeop
    else:
        ZorH = (RE * Hgeop) / (RE + Hgeop)

    # Unit Conversion OUT
    rho = Q_(rho, ureg.kg/ureg.m**3)
    a = Q_(a, ureg.m/ureg.s)
    temp = Q_(temp, ureg.degK)
    press = Q_(press, ureg.Pa)
    kvisc = Q_(kvisc, ureg.m**2/ureg.s)
    ZorH = Q_(ZorH, ureg.m)
    if Units[1] == 'US':
        rho = rho.to(ureg.slug/ureg.ft**3)
        a = a.to(ureg.ft/ureg.s)
        temp = temp.to(ureg.degR)
        press = press.to(ureg.psi)
        kvisc = kvisc.to(ureg.ft**2/ureg.s)
        ZorH = ZorH.to(ureg.ft)
    if not dimensioned:
        rho = rho.magnitude
        a = a.magnitude
        temp = temp.magnitude
        press = press.magnitude
        kvisc = kvisc.magnitude
        ZorH = ZorH.magnitude
    return rho, a, temp, press, kvisc, ZorH


# Example usage
if __name__ == "__main__":
    Rho, A, Temp, Press, Kvisc, Z_orH = atmos(
        H_in=0, T_offset=0, Units='SI', GeomFlag=False)

    print("Density:", f"{Rho}")
    print("Speed of sound:", f"{A}")
    print("Temperature:", f"{Temp}")
    print("Pressure:", f"{Press}")
    print("Kinematic viscosity:", f"{Kvisc}")
    print("Height or altitude:", f"{Z_orH}")
