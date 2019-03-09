import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col


# Elementary transfer functions

def make_V_func(V):
    '''
    Returns the transfer function F(s) = V.
    '''
    def V_func(s):
        return V
    return V_func

def make_D_func(T_D):
    '''
    Returns the transfer function F(s) = T_d * s.
    '''
    def D_func(s):
        return T_D * s
    return D_func

def make_I_func(T_I):
    '''
    Returns the transfer function F(s) = 1 / (T_I * s).
    '''
    def I_func(s):
        return 1.0 / (T_I * s)
    return I_func

def make_PT1_func(T):
    '''
    Returns the transfer function F(s) = 1 / (T*s + 1).
    '''
    def PT1_func(s):
        return 1.0 / (T * s + 1)
    return PT1_func

def make_PT2_func(omega, D):
    '''
    Returns the transfer function F(s) = 1 / (1/omega^2 * s^2 + 2D/omega * s + 1).
    '''
    def PT2_func(s):
        return 1.0 / ((s * s) / (omega * omega) + (2 * D * s) / omega + 1)
    return PT2_func

def make_DT1_func(T):
    '''
    Returns the transfer function F(s) = T*s + 1.
    '''
    def DT1_func(s):
        return T * s + 1
    return DT1_func


# Composite transfer functions

def make_prod_func(functions):
    '''
    Returns the product of the given transfer functions.
    '''
    def prod_func(s):
        prod = 1
        for F in functions:
            prod *= F(s)
        return prod
    return prod_func

def make_sum_func(functions):
    '''
    Returns the sum of the given transfer functions.
    '''
    def sum_func(s):
        sum = 0
        for F in functions:
            sum += F(s)
        return sum
    return sum_func

def make_feedback_func(F_forward, F_0):
    '''
    Returns the transfer function of the feedback loop F(s) = F_forward / (1 + F_0)
    '''
    def feedback_func(s):
        return F_forward(s) / (1 + F_0(s))
    return feedback_func


# Linear, logarithmic and hue intervalls

def intervall(a=0.0, b=1.0, delta=0.01):
    '''
    Returns the np.linspace from a to b with step delta.
    '''
    return np.linspace(a, b, (b - a) / delta + 1)

def log_intervall(a=1.0, b=10.0, num=50.0, base=10.0):
    '''
    Returns the np.logspace from base**a to base**b with num samples per exponent increase.
    '''
    return np.logspace(a, b, num=(b-a)*num, base=base)

def hue_intervall(num, saturation, value, start=0.0):
    '''
    Returns a list of rgb triples visually evenly spaced in regard of hues; start markes first color hue in degrees.
    '''
    start /= 360.0

    hues = np.arange(0.0, 1.0, 1.0 / num)
    for i in range(len(hues)):
        hues[i] = math.pow(hues[i], 1.5)
    
    colors = []
    for hue in hues:
        hsv = ((hue + start) % 1.0, saturation, value)
        colors.append(col.hsv_to_rgb(hsv))
       
    return colors
    

# Bode diagramm

def make_dB_func(function):
    '''
    Returns the absolute value in dB of the given transfer function.
    '''
    def dB_func(omega):
        return 20 * np.log10(abs(function(1j * omega)))
    return dB_func

def make_phase_func(function, phi_min, phi_max):
    '''
    Returns the phase value in degrees in the range between phi_min and phi_max of the given transfer function.
    '''
    def phase_func(omega):
        phi = math.degrees(cmath.phase(function(1j * omega)))
        if phi > phi_max:
            phi -= 360
        if phi < phi_min:
            phi += 360
        return phi
    return phase_func

class EmpytBode(object):
    '''
    Empyt bode diagramm for hand sketches.
    '''
    def __init__(self, start_exp, end_exp, db_min, db_max, db_delta, phi_min, phi_max, phi_delta):
        '''
        Creates intervalls for frequency omega, absolute value in dB and phase phi.
        '''
        self.omega = log_intervall(a=start_exp, b=end_exp)
        self.db_lim = intervall(db_min, db_max, db_delta)
        self.phi_lim = intervall(phi_min, phi_max, phi_delta)

    def plot(self):
        '''
        Returns matplotlib figure of the bode diagramm.
        '''
        fig, ax_db = plt.subplots(1, 1)
        fig.subplots_adjust(left=0.125, right=0.875, bottom=0.125)
        ax_phi = ax_db.twinx()

        ax_db.set_xscale("log")
        ax_phi.set_xscale("log")

        ax_db.set_ylim([self.db_lim[0], self.db_lim[-1]])
        ax_db.set_yticks(self.db_lim)

        ax_db.set_xlim(self.omega[0], self.omega[-1])

        ax_phi.set_ylim([self.phi_lim[0], self.phi_lim[-1]])
        ax_phi.set_yticks(self.phi_lim)

        ax_db.set_xlabel(r"$\omega \ / \ \frac{1}{s}$")
        ax_db.set_ylabel("Betrag / dB")
        ax_phi.set_ylabel("Phase / °")

        ax_db.grid(b=True, which="both", axis="both")

        return fig

    def save(self, path="", filename="plot.png"):
        '''
        Creates and saves bode diagramm at path/filename.
        '''
        self.plot().savefig(path + filename)

    def show(self):
        '''
        Creates and shows bode diagramm.
        '''
        self.plot().show()

class Bode(object):
    '''
    Bode diagramm of an arbitrary number of transfer functions.
    '''
    def __init__(self, functions, labels, start_exp, end_exp, db_min, db_max, db_delta, phi_min, phi_max, phi_delta):
        '''
        Creates intervalls for frequency omega, absolute values in dB and phases phi.
        Sets functions, labels and colors as attributes.
        '''
        self.functions = functions
        self.labels = labels
        self.colors = hue_intervall(len(functions), 1, 0.8)

        self.omega = log_intervall(a=start_exp, b=end_exp)
        self.dB = []
        self.phi = []

        for F in functions:
            dB_func = make_dB_func(F)
            self.dB.append([dB_func(f) for f in self.omega])
            
            phase_func = make_phase_func(F, phi_min, phi_max)
            self.phi.append([phase_func(f) for f in self.omega])
        
        self.db_lim = intervall(db_min, db_max, db_delta)
        self.phi_lim = intervall(phi_min, phi_max, phi_delta)

    def plot(self):
        '''
        Returns matplotlib figure of the bode diagramm.
        '''
        fig, ax_db = plt.subplots(1, 1)
        fig.subplots_adjust(left=0.125, right=0.875, bottom=0.125)
        ax_phi = ax_db.twinx()

        ax_db.set_xscale("log")
        ax_phi.set_xscale("log")

        for i in range(0, len(self.dB)):
            ax_db.plot(self.omega, self.dB[i], label=r"$|$" + self.labels[i] + r"$|$", color=self.colors[i], linewidth=2)
            ax_phi.plot(self.omega, self.phi[i], label=r"$\varphi ($" + self.labels[i] + r"$)$", color=self.colors[i], linewidth=2, linestyle=":")

        ax_db.set_ylim([self.db_lim[0], self.db_lim[-1]])
        ax_db.set_yticks(self.db_lim)

        ax_db.set_xlim(self.omega[0], self.omega[-1])

        ax_phi.set_ylim([self.phi_lim[0], self.phi_lim[-1]])
        ax_phi.set_yticks(self.phi_lim)

        #ax.set_xlim(xlim)
        #ax.set_ylim(ylim)

        ax_db.set_xlabel(r"$\omega \ / \ \frac{1}{s}$")
        ax_db.set_ylabel("Betrag / dB")
        ax_phi.set_ylabel("Phase / °")

        #ax.set_title(title)

        ax_db.legend(loc="upper left")
        ax_phi.legend(loc="upper right")

        ax_db.grid(b=True, which="both", axis="both")
        #ax_phi.grid(b=True, which="major", axis="y")

        return fig

    def save(self, path="", filename="plot.png"):
        '''
        Creates and saves bode diagramm at path/filename.
        '''
        self.plot().savefig(path + filename)

    def show(self):
        '''
        Creates and shows bode diagramm.
        '''
        self.plot().show()
