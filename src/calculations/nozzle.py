import calculations.init as init
from scipy.optimize import fsolve

# # ------------ Constants ------------------
gamma = init.GAMMA
# # ------------------------------------------

def area_mach(mach):
        """
        Calculate the area ratio for a given Mach number using the isentropic flow relations.
        """
        area_ratio = (1 / mach) * ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * mach**2))**((gamma + 1) / (2 * (gamma - 1)))
        return area_ratio


def mach_design_exit(area_ratio):
    """
    Calculate the design Mach number at the nozzle exit for a given area ratio.
    """
    function = lambda mach_exit: area_mach(mach_exit) - area_ratio
    mach_exit = fsolve(function, 2.0)[0]
    return mach_exit


def static_pressure_ratio(mach):
        """
        Calculate the static pressure ratio for a given Mach number.
        """
        p_p0 = (1 + ((gamma - 1) / 2) * mach**2) ** (-gamma / (gamma - 1))
        return p_p0