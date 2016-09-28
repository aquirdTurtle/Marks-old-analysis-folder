# -*- coding: utf-8 -*-
"""
Created on Fri May  6 16:17:37 2016
This script is for the calculation of the trapping frequency for a 1D Optical lattice. It could be expanded if someone wished...
@author: Mark
"""
import math
""""""
""" Input Constants """
""""""
# in kg
mass = 1.443*10**(-25);
# in meters
trapLightWavelength = 820*10**(-9);
#trapLightWavelength = 915*10**(-9);
# in hertz
#desiredTrapFrequency = 2*math.pi*75*10**3;
# from Markus's thesis
desiredTrapFrequency = 2*math.pi*100*10**3;
# From other thesis
#desiredTrapFrequency = 2*math.pi*44*10**3;
# in meters
beamWaistAtAtoms = 170*10**(-6);
#beamWaistAtAtoms = 120*10**(-6);
""""""
""" Things that I don't change """
""""""
# in meters/second
speedOfLight = 2.998*10**8
# This is different from the wavelength because of the bowtie configuration.
# in meters
lightWavelengthAtTrap = math.sqrt(2)*trapLightWavelength;
#lightWavelengthAtTrap = trapLightWavelength;
# in Hz
d1TransitionFrequency = 384.2*10**12;
d2TransitionFrequency = 377.1*10**12;
d1LineWidth = 36.1*10**6;
d2LineWidth = 38.11*10**6;
detuningFromD1 = speedOfLight / trapLightWavelength - d1TransitionFrequency;
detuningFromD2 = speedOfLight / trapLightWavelength - d2TransitionFrequency;
power = -(((2*math.pi)**2 * mass * (beamWaistAtAtoms**2) * (desiredTrapFrequency**2)*(lightWavelengthAtTrap**2))
        / (8 * (speedOfLight**2)* (2 * d2LineWidth/(d2TransitionFrequency**3*detuningFromD2) + (d1LineWidth/(d1TransitionFrequency**3*detuningFromD1)))));
print(power);

""""""
""" TWEEZER TRAP """
""""""
"""
Power = 0.2
Depth = -((math.pi*speedOfLight**2)*averageNaturalLineWidth*((2/detuningFromD2) + (1/detuningFromD1)) / (2*averageTransitionFrequency**3))*2*Power / (math.pi*beamWaistAtAtoms**2)
print(Depth / ((1.054*10**(-34))**2*((2*math.pi / trapLightWavelength))**2/(2*mass)))
frequency = math.sqrt(4*Depth / (mass * beamWaistAtAtoms**2))
print(frequency)
"""