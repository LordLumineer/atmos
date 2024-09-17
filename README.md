# atmos
ATMOS Find gas properties in earth's atmosphere.

This module provides a function `atmos` which computes the gas properties in
earth's atmosphere. The function takes as input the geopotential altitude in
meters and a temperature offset in degrees Celsius. It can also take as option
the unit system for inputs and outputs and the use of geometric instead of
geopotential altitude.

The function returns a tuple containing the density, speed of sound,
temperature, pressure, kinematic viscosity and height or altitude of the
atmosphere at the given input.

## Parameters

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

## Returns

- rho : float or pint.Quantity
  - Density in kilograms per cubic meter or slugs per cubic foot.
- a : float or pint.Quantity
  - Speed of sound in meters per second or feet per second.
- T : float or pint.Quantity
  - Temperature in degrees Kelvin or Rankine.
- P : float or pint.Quantity
  - Pressure in pascals or pounds per square foot.
- kvisc : float or pint.Quantity
  - Kinematic viscosity in square meters per second or square feet per second.
- ZorH : float or pint.Quantity
  - Height or altitude in meters or feet.

## Example

```python
>>> from atmos import atmos
>>> rho, a, T, P, kvisc, ZorH = atmos(H_in=1000, T_offset=0, Units='SI', GeomFlag=False)
>>> print("Density:", rho)
>>> print("Speed of sound:", a)
>>> print("Temperature:", T)
>>> print("Pressure:", P)
>>> print("Kinematic viscosity:", kvisc)
>>> print("Height or altitude:", ZorH)
```

## Notes

ATMOS Find gas properties in earth's atmosphere.

ATMOS by itself gives the atmospheric properties at sea level on a standard day.

ATMOS(H) returns the properties of the 1976 Standard Atmosphere at geopotential altitude H (meters).

ATMOS(H,dT) returns properties when the temperature is dT degrees offset from standard conditions.

ATMOS(H,dT,Units) specifies units for the inputs outputs.

- The unit system can be specified as a string or a list of strings. 
- If a single string is provided, it will be used for both the input and output units.  
- If a list of two strings is provided, the first string will be used for the input units and the second string be used for the output units.  
- The supported unit systems are: SI (default), US (a.k.a. Imperial, English), (others pint supported unit systems are not yet supported).  
- If a pint.Quantity is provided as an input, the units will be automatically detected.  

For example:

- To use the SI system, set Units to '', or 'SI'.
- To use the US system, set Units to 'US'.
- To use different units for the input and output, pass a list of strings, e.g. ['US', 'SI'].

The temperature offset (dT) is in the same units as the input temperature, so when converting between Celsius and Fahrenheit, use only the scaling factor (dC/dF = dK/dR = 5/9).

Units are as follows:

| Input                   | SI (default) | US        |
| ----------------------- | ------------ | --------- |
| H:  Altitude            | m            | ft        |
| dT: Temp. offset        | K / Δ°C      | °R / Δ°F  |

| Output                  | SI (default) | US        |
| ----------------------- | ------------ | --------- |
| rho: Density            | kg/m^3       | slug/ft^3 |
| a: Speed of sound       | m/s          | ft/s      |
| T: Temperature          | K            | °R        |
| P: Pressure             | Pa           | psi       |
| nu: Kinem. viscosity    | m^2/s        | ft^2/s    |
| orH: Height or altitude | m            | ft        |

This atmospheric model is not recommended for use at altitudes above 86 km geometric height (84852 m/278386 ft geopotential) and returns NaN for altitudes above 100 km geopotential.

References: ESDU 77022; [www.pdas.com/atmos.html](www.pdas.com/atmos.html)
