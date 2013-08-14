#!/uar/bin/ewnv python

"""
OilWeathering.py

module to compute "pseudo component" weathering

primarily the weather_curve class

Built-in Oil Types are in the OilTypes dict.

NOTE:

To compute the half lives of components subject to multiple processes (linear),
such as evaporation and bio-degradation, you can use the following formula:

H_t = 1 / ( 1/t_w + 1/t_b)


"""

import numpy as np

class weather_curve:
    def __init__(self, C, H):
        """
        Simple weathering computation for a three "component" oil
        
        Each component is a fraction of the total mass and has its own half-life
        
        (C1, C2, C3, ... Ci) are the fractions of each component (must add up to 1.0)
        (H1, H2, H3, ....Hi) are the half lives of each component (in hours)
        
        """
        self.C = np.asarray(C, dtype=np.float32).reshape(-1,)
        self.H = np.asarray(H, dtype=np.float32).reshape(-1,)
        if round(self.C.sum(), 6)  != 1.0: # only six digit, because float32
            raise ValueError("The three constants must add up to one. These add up to: %f"%self.C.sum())
        if len(self.H) != len(self.C):
            raise ValueError("There must be the same number of component fractions as half lives")
            
    
    def weather (self, M_0, time):
        """
        compute what mass is left at time specified
        
        M_0 is the initial mass
        time is the time from release, in hours
        
        returns the mass remaining
        
        """
        M_0 = np.asarray(M_0, dtype=np.float32)
        time = np.asarray(time, dtype=np.float32)
        half = np.float32(0.5)
        
        M = 0
        for C, H in zip(self.C, self.H):
            M += C * half**(time/H)
        return M_0 * M


## Parameters for combined weathering and bio-degradation for "medium crude"
## used for FL Staits TAP analysis
mass_fractions =       [0.25,   0.1, 0.107,    0.2,  0.186,   0.109,    0.048]
combined_half_lives =  [21.0, 422.0,   2.0, 1358.0, 1982.0,  7198.0,  14391.0]
    
OilTypes = {None: None,
            # Medium Crude parameters from OSSM
            'MediumCrude': weather_curve( ( .22,  .26, .52),
                                          (14.4, 48.6, 1e9)),
            "FL_Straits_MediumCrude": weather_curve(mass_fractions, combined_half_lives),
            }

