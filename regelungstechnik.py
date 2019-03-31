#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col


# Matplotlib settings

plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["text.usetex"] = True


# Metric prefixes

def unify(val):
    prefix = ["",
              "k", "M", "G", "T", "P", "E", "Z", "Y",
              "y", "z", "a", "f", "p", "n", "µ", "m"]

    i = 0

    val = np.asarray(val, dtype=np.float64)

    # Multiple
    if val.max() > 1:
        while val.max() >= 1000:
            val *= 1e-3
            i += 1

    # Fraction
    elif val.max() < 1:
        while val.max() < 1:
            val *= 1e3
            i -= 1

    if val.size == 1:
        val = float(val)

    return val, prefix[i]


# String representations

def str_V(V, lin=True, dB=True, prefix=False, text="V"):
    """
    Convert V into a string representation.
    Choose if you want the linear or dB value.
    Choose if V shall get its metric prefix.
    The displayed name can be specified in text.
    """
    if prefix:
        V, pre = unify(V)
    else:
        pre = ""

    val = text

    if lin:
        val += " = " + str(V) + pre

    if dB:
        val += " = " + str(dB(V)) + "dB"

    return val

def str_T(T, prefix=True, text="T"):
    """
    Convert T into a string representation.
    Choose if T shall get its metric prefix.
    The displayed name can be specified in text.
    """
    if prefix:
        T, pre = unify(T)
    else:
        pre = ""

    val = text
    val += " = " + str(T) + pre + "s"

    return val

def str_omega(omega, prefix=True, text="omega"):
    """
    Convert omega into a string representation.
    Choose if omega shall get its metric prefix.
    The displayed name can be specified in text.
    """
    if prefix:
        omega, pre = unify(omega)
    else:
        pre = ""

    val = text
    val += " = " + str(omega) + pre + "/s"

    return val

def str_f(f, prefix=True, text="f"):
    """
    Convert f into a string representation.
    Choose if f shall get its metric prefix.
    The displayed name can be specified in text.
    """
    if prefix:
        omega, pre = unify(omega)
    else:
        pre = ""

    val = text
    val += " = " + str(omega) + pre + "Hz"

    return val


# Basic functions and transformations

def delta(t):
    val = np.zeros(t.size, dtype=np.float64)
    val[t == 0] = 1
    return val

def step(t):
    val = np.ones(t.size, dtype=np.float64)
    val[t < 0] = 0
    return val

def ramp(t):
    val = np.copy(t)
    val[t < 0] = 0
    return val


# Fourier transform

def fft(y):
    l = y.size // N + 1
    val = np.fft.fft(y)
    val[l:] = 0
    return val

def ifft(Y):
    l = Y.size // N + 1
    val = np.fft.irfft(Y[:l])
    return val


# Basic elements

class BasicElement(object):
    def __init__(self):
        """
        Template for a basic element.
        """
        def H(s):
            """
            The transfer function H(s).
            """
            s = np.asarray(s, dtype=np.complex64)
            val = None
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t).
            """
            t = np.asarray(t, dtype=np.float64)
            val = None
            return val
        self.h = h

        def w(t):
            """
            The step response w(t).
            """
            t = np.asarray(t, dtype=np.float64)
            val = None
            return val
        self.w = w

        self.roots = None
        self.poles = None

        self.counter = None
        self.denominator = None

    def __str__(self):
        val = "Basic Element"
        return val


class P(BasicElement):
    def __init__(self, V=1, dB=False):
        """
        The basic element P (proportional).
        V may be given in dB.
        """
        if dB:
            V = lin(V)

        self.V = V

        def H(s):
            """
            The transfer function H(s) = V.
            """
            s = np.asarray(s, dtype=np.complex64)
            val = V * np.ones(shape=s.shape, dtype=np.complex64)
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t) = [V, 0, ...]
            """
            t = np.asarray(t, dtype=np.float64)
            val = np.zeros(shape=t.shape, dtype=np.float64)
            val[0] = V
            return val
        self.h = h

        def w(t):
            """
            The step response w(t) = V
            """
            t = np.asarray(t, dtype=np.float64)
            val = V * np.ones(shape=t.shape, dtype=np.float64)
            return val
        self.w = w

        self.roots = []
        self.poles = []

        self.counter = V
        self.denominator = 1

    def __str__(self):
        val = "Gain P"
        val += " with " + str_V(self.V)
        return val


class I(BasicElement):
    def __init__(self, T=1):
        """
        The basic element I (integrator).
        """
        self.T = T

        def H(s):
            """
            The transfer function H(s) = 1 / Ts.
            """
            s = np.asarray(s, dtype=np.complex64)
            val = 1 / (T * s)
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t) = V.
            """
            t = np.asarray(t, dtype=np.float64)
            val = V * np.ones(shape=t.shape, dtype=np.float64)
            return val
        self.h = h

        def w(t):
            """
            The step response w(t) = 1 / T * t.
            """
            t = np.asarray(t, dtype=np.float64)
            val = 1 / T * t
            return val
        self.w = w

        self.roots = []
        self.poles = [0 + 0j]

        self.counter = 1
        self.denominator = D(T).H

    def __str__(self):
        val = "Integrator I"
        val += " with " + str_T(self.T)
        return val


class D(BasicElement):
    def __init__(self, T=1):
        """
        The basic element D (differentiator).
        """
        self.T = T

        def H(s):
            """
            The transfer function H(s) = T * s.
            """
            s = np.asarray(s, dtype=np.complex64)
            val = T * s
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t) = 0.
            """
            t = np.asarray(t, dtype=np.float64)
            val = np.zeros(shape=t.shape, dtype=np.float64)
            return val
        self.h = h

        def w(t):
            """
            The step response w(t) = [inf, 0, ...].
            """
            t = np.asarray(t, dtype=np.float64)
            val = np.zeros(shape=t.shape, dtype=np.float64)
            val[0] = np.inf
            return val
        self.w = w

        self.roots = [0 + 0j]
        self.poles = []

        self.counter = H
        self.denominator = 1

    def __str__(self):
        val = "Differentiator D"
        val += " with " + str_T(self.T)
        return val


class PT1(BasicElement):
    def __init__(self, T=1, V=1, dB=False):
        """
        The basic element PT1 (low pass of order 1).
        V may be given in dB.
        """
        if dB:
            V = lin(V)

        self.T = T
        self.V = V

        def H(s):
            """
            The transfer function H(s) = V / (Ts + 1).
            """
            s = np.asarray(s, dtype=np.complex64)
            val = V / (T * s + 1)
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t).
            """
            t = np.asarray(t, dtype=np.float64)
            val = V / T * np.exp(-1 * t / T)
            return val
        self.h = h

        def w(t):
            """
            The step response w(t).
            """
            t = np.asarray(t, dtype=np.float64)
            val = V * (1 - np.exp(-1 * t / T))
            return val
        self.w = w

        self.roots = []
        self.poles = [-1 / T]

        self.counter = V
        self.denominator = PD1(T=T, V=V, dB=dB).H

    def __str__(self):
        val = "Low pass PT1 of order 1"
        val += " with " + str_V(self.V)
        val += " and " + str_T(self.T)
        return val


class PT2(BasicElement):
    def __init__(self, omega=1, D=1, V=1, dB=False):
        """
        The basic element PT2 (low pass of order 2).
        V may be given in dB.
        """
        if D < 0:
            print("Error: Attenuation must not be negative.")
            return

        if dB:
            V = lin(V)

        self.omega = omega
        self.D = D
        self.V = V

        def H(s):
            """
            The transfer function H(s) = V / (s/omega ** 2 + 2D/omega * s + 1).
            """
            s = np.asarray(s, dtype=np.complex64)
            val = V / (s / omega ** 2 + 2 * D / omega * s + 1)
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t) of a PT2 at different attenuations.
            """
            t = np.asarray(t, dtype=np.float64)
            if D == 0:
                val = V / (2 * omega) * np.sin(omega * t)
            if D > 0 and D < 1:
                val = V / (2 * omega * np.sqrt(1 - D*D))
                val *= np.exp(-1 * D * omega * t)
                val *= np.sin(np.sqrt(1 - D*D) * omega * t)
            elif D == 1:
                val = V * omega * (2 - omega * t) * np.exp(-1 * omega * t)
            else:
                D1 = np.sqrt(D ** 2 - 1)
                val = -1 * V * omega / (2 * D1)
                val *= np.exp(-1 * (D - D1) * omega * t)
            val[t < 0] = 0
            return val
        self.h = h

        def w(t):
            """
            The step response w(t) of a PT2 at different attenuations.
            """
            t = np.asarray(t, dtype=np.float64)
            if D == 0:
                val = V - V * np.cos(omega * t)
            elif D > 0 and D < 1:
                theta = np.arccos(D)
                D1 = np.sqrt(1 - D ** 2)
                val = -1 * V * np.exp(-1 * omega * t) / D1
                val *= np.sin(D1 * omega * t + theta)
                val += V
            elif D == 1:
                val = V - V * (1 - omega * t) * np.exp(-1 * omega * t)
            else:
                D1 = np.sqrt(D ** 2 - 1)
                val = V / (2 * (D - D1) * D1)
                val *= np.exp(-1 * (D - D1) * omega * t)
            val[t < 0] = 0
            return val
        self.w = w

        self.roots = []

        if D == 0:
            self.poles = [0 + 1j * omega, 0 - 1j * omega]
        elif D > 0 and D < 1:
            real = -1 * omega * D
            imag = 1j * omega * np.sqrt(1 - D ** 2)
            self.poles = [real + imag, real - imag]
        elif D == 1:
            self.poles = [-1 * omega]
        else:
            p1 = -1 * omega * D
            p2 = omega * np.sqrt(D ** 1 - 1)
            self.poles = [p1 - p2, p1 + p2]

        self.counter = V
        self.denominator = PD2(omega=omega, D=D, V=V, dB=dB).H

    def __str__(self):
        val = "Low pass PT2 of order 2"
        val += " with " + str_V(self.V)
        val += " and " + str_omega(self.omega)
        val += " and " + str_V(self.D, dB=False, text="D")
        return val


class PD1(BasicElement):
    def __init__(T=1, V=1, dB=False):
        """
        The basic element PD1 (allowance of order 1).
        V may be given in dB.
        """
        if dB:
            V = lin(V)

        self.T = T
        self.V = V

        def H(s):
            """
            The transfer function H(s) = V * (Ts + 1).
            """
            s = np.asarray(s, dtype=np.complex64)
            val = V * (T * s + 1)
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t).
            """
            t = np.asarray(t, dtype=np.float64)
            val = V * delta(t)
            return val
        self.h = h

        def w(t):
            """
            The step response w(t).
            """
            t = np.asarray(t, dtype=np.float64)
            val = V * (T * delta(t) + step(t))
            return val
        self.w = w

        self.roots = [-1 / T]
        self.poles = []

        self.counter = H
        self.denominator = 1

    def __str__(self):
        val = "Allowance PD1 of order 1"
        val += " with " + str_V(self.V)
        val += " and " + str_T(self.T)
        return val


class PD2(BasicElement):
    def __init__(self, omega=1, D=1, V=1, dB=False):
        """
        The basic element PT2 (allowance of order 2).
        V may be given in dB.
        """
        if D < 0:
            print("Error: Attenuation must not be negative.")
            return

        if dB:
            V = lin(V)

        self.omega = omega
        self.D = D
        self.V = V

        def H(s):
            """
            The transfer function H(s) = V / (s/omega ** 2 + 2D/omega * s + 1).
            """
            s = np.asarray(s, dtype=np.complex64)
            val = V * (s / omega ** 2 + 2 * D / omega * s + 1)
            return val
        self.H = H

        def h(t):
            """
            The impulse response h(t) of a PT2 at different attenuations.
            """
            t = np.asarray(t, dtype=np.float64)
            val = V * delta(t)
            return val
        self.h = h

        def w(t):
            """
            The step response w(t) of a PT2 at different attenuations.
            """
            t = np.asarray(t, dtype=np.float64)
            val = V * (2 * D / omega * delta(t) + step(t))
            return val
        self.w = w


        if D == 0:
            self.roots = [0 + 1j * omega, 0 - 1j * omega]
        elif D > 0 and D < 1:
            real = -1 * omega * D
            imag = 1j * omega * np.sqrt(1 - D ** 2)
            self.roots = [real + imag, real - imag]
        elif D == 1:
            self.roots = [-1 * omega]
        else:
            p1 = -1 * omega * D
            p2 = omega * np.sqrt(D ** 1 - 1)
            self.roots = [p1 - p2, p1 + p2]

        self.poles = []

        self.counter = H
        self.denominator = 1

    def __str__(self):
        val = "Allowance PD2 of order 2"
        val += " with " + str_V(self.V)
        val += " and " + str_omega(self.omega)
        val += " and " + str_V(self.D, dB=False, text="D")
        return val


# Composite transfer functions

def PROD(functions):
    """
    Returns the product of the given transfer functions
    as a new transfer function.
    """
    def F(s, g=[]):
        s = np.asarray(s, dtype=np.complex64)
        g = np.asarray(g, dtype=np.float64)
        prod = np.ones(s.size, dtype=np.complex64)
        imp_res = np.zeros(s.size, dtype=np.float64)
        imp_res[0] = 1
        for F in functions:
            prod *= F(s, g=g)

        return prod
    return F


def SUM(functions):
    """
    Returns the sum of the given transfer functions
    as a new transfer function.
    """
    def F(s, g=[]):
        s = np.asarray(s, dtype=np.complex64)
        sum = np.zeros(s.shape, dtype=np.complex64)
        for F in functions:
            sum += F(s, g=[])
        return sum
    return F


def FEEDBACK_LOOP(F1, F2):
    """
    Returns the transfer function of the feedback loop
    F(s, g=[]) = F1 / (1 + F1 * F2)
    """
    def F(s, g=[]):
        s = np.asarray(s, dtype=np.complex64)
        return F1(s) / (1 + F1(s) * F2(s))
    return F


# Evaluate abs and phase of complex values in dB and degrees

def dB(val):
    """
    Returns the absolute value in dB of the given values.
    """
    return 20 * np.log10(np.abs(val))


def lin(val):
    """
    Returns the linear values of the given dB values.
    """
    return 10 ** (val / 20)


def phase(val, min, max):
    """
    Returns the phase value in degrees in the range between
    v_min and v_max of the given values.
    """
    phi = np.angle(val, deg=True)
    # Move phi down
    while (phi > max).any():
        phi[phi > max] -= 360
    # Move phi up
    while (phi < min).any():
        phi[phi < min] += 360
    # Move phi down as far as possible
    while (phi - 360 > min).any():
        phi[phi - 360 > min] -= 360
    return phi


# Hue color intervall

def hue_intervall(num, saturation, value, start=0.0):
    """
    Returns a list of rgb triples visually evenly spaced in regard of hues;
    start markes first color hue in degrees.
    """
    # Matplotlib uses hues in range [0,1].
    start /= 360.0

    hues = np.arange(0.0, 1.0, 1.0 / num)
    # The hues in the hsv color space are visually not evenly distributed.
    # To compensate this effect, we calculate hue**1.5.
    for i in range(len(hues)):
        hues[i] = math.pow(hues[i], 1.5)

    colors = []
    for hue in hues:
        hsv = ((hue + start) % 1.0, saturation, value)
        colors.append(col.hsv_to_rgb(hsv))

    return colors


# Bode diagramm

class BodeDiagram(object):
    """
    Bode diagramm of an arbitrary number of transfer functions.
    """
    def __init__(self, functions, labels, start, stop, ticks,
        dB_delta=20, phi_delta=45, N=1024):
        """
        Takes a list of transfer functions with corresponding labels and
        creates a bode diagramm from 10**start to 10**stop and a given list
        of ticks.
        """
        self.functions = functions
        self.labels = labels
        self.colors = hue_intervall(len(functions), 1, 0.8)

        self.omega = np.logspace(start=start, stop=stop, num=N)

        self.dBs = []
        self.phis = []

        self.db_ticks = dB_delta * np.asarray(ticks)
        self.phi_ticks = phi_delta * np.asarray(ticks)

        for F in functions:
            values = F(1j * self.omega)
            phi_min, phi_max = self.phi_ticks[0], self.phi_ticks[-1]
            self.dBs.append(dB(values))
            self.phis.append(phase(values, phi_min, phi_max))


    def plot(self, pick=None):
        """
        Returns matplotlib figure of the bode diagramm.
        """
        fig, ax_db = plt.subplots(figsize=(8, 4.5), dpi=240)
        fig.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.95)
        ax_phi = ax_db.twinx()

        ax_db.set_xscale("log")
        ax_phi.set_xscale("log")

        dBs = np.asarray(self.dBs)
        phis = np.asarray(self.phis)
        labels = np.asarray(self.labels)
        colors = np.asarray(self.colors)

        canvas = False

        if pick is not None:
            pick = np.asarray(pick, dtype=int)
            if pick.size == 0:
                canvas = True
            dBs = dBs[pick]
            phis = phis[pick]
            labels = labels[pick]
            colors = colors[pick]

        for dB, phi, label, color in zip(dBs, phis, labels, colors):
            dB_label = r"$|$" + label + r"$|$"
            phi_label = r"$\varphi ($" + label + r"$)$"

            ax_db.plot(self.omega, dB, label=dB_label,
                       color=color, linewidth=2)
            ax_phi.plot(self.omega, phi, label=phi_label,
                        color=color, linewidth=2, linestyle=":")

        ax_db.set_xlim(self.omega[0], self.omega[-1])

        ax_db.set_ylim([self.db_ticks[0], self.db_ticks[-1]])
        ax_db.set_yticks(self.db_ticks)

        ax_phi.set_ylim([self.phi_ticks[0], self.phi_ticks[-1]])
        ax_phi.set_yticks(self.phi_ticks)

        ax_db.set_xlabel(r"$\omega \ / \ \frac{1}{s}$")
        ax_db.set_ylabel("Betrag / dB")
        ax_phi.set_ylabel("Phase / °")

        if not canvas:
            ax_db.legend(loc="upper left")
            ax_phi.legend(loc="upper right")

        ax_db.grid(b=True, which="both", axis="both")

        return fig


    def save(self, pick=None, path="", filename="plot.png"):
        """
        Creates and saves bode diagramm at path/filename.
        """
        fig = self.plot(pick=pick)
        fig.savefig(path + filename)


    def show(self, pick=None):
        """
        Creates and shows bode diagramm.
        """
        fig = self.plot(pick=pick)
        fig.show()


# Impulse response

class StepResponse(object):
    """
    Step response of an arbitrary number of transfer functions.
    """
    def __init__(self, functions, labels, duration, N=1024):
        """
        Sets functions, labels and colors as attributes.
        Computes sample rate, time and omega from duration and N.
        """
        self.functions = functions
        self.labels = labels
        self.colors = hue_intervall(len(functions), 1, 0.8)

        sample_rate = N * 2 * np.pi / duration
        self.time = np.linspace(0, duration, N)
        omega = np.linspace(0, sample_rate, N)

        self.steps = []

        for F in functions:
            # Inverse real fourier g.shape == s.shape, only takes lower half of spectrum
            spectrum = F(1j * omega)
            spectrum = spectrum[0 : N // 2 + 1]
            impulse_response = np.fft.irfft(spectrum)

            # Step response = Integrated impulse response
            step_response = impulse_response.cumsum()

            self.steps.append(step_response)


    def plot(self, pick=None, v_min=None, v_max=None):
        """
        Returns matplotlib figure of the step responses.
        """
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=240)
        fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.95)

        steps = np.asarray(self.steps)
        labels = np.asarray(self.labels)
        colors = np.asarray(self.colors)

        canvas = False

        if pick is not None:
            pick = np.asarray(pick, dtype=int)
            if pick.size == 0:
                canvas = True
            steps = steps[pick]
            labels = labels[pick]
            colors = colors[pick]

        for step, label, color in zip(steps, labels, colors):
            ax.plot(self.time, step, label=label, color=color, linewidth=2)

        ax.set_xlim(self.time[0], self.time[-1])

        if v_min is None:
            v_min = 0 if canvas else min([step.min() for step in steps])
        if v_max is None:
            v_max = 1 if canvas else max([step.max() for step in steps])

        ax.set_ylim(v_min, v_max)

        ax.set_xlabel(r"$t \ / \ s$")
        ax.set_ylabel("Regelgröße")

        if not canvas:
            ax.legend(loc="upper right")

        ax.grid(b=True, which="both", axis="both")

        return fig


    def save(self, pick=None, v_min=None, v_max=None,
             path="", filename="plot.png"):
        """
        Creates and saves step response at path/filename.
        """
        fig = self.plot(pick=pick, v_min=v_min, v_max=v_max)
        fig.savefig(path + filename)


    def show(self, pick=None, v_min=None, v_max=None):
        """
        Creates and shows step response.
        """
        self.plot(pick=pick, v_min=v_min, v_max=v_max).show()
