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


# Elementary transfer functions

def P(V, dB=False):
    """
    Returns the transfer function F(s) = V.
    V may be given in dB.
    """
    if dB:
        V = lin(V)
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        return V * np.ones(shape=s.shape, dtype=np.complex64)
    return F


def D(T):
    """
    Returns the transfer function F(s) = T * s.
    """
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        return T * s
    return F


def I(T):
    """
    Returns the transfer function F(s) = 1 / Ts.
    """
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        duration = 2j * np.pi / (s[1] - s[0])
        val = duration / T * np.ones(s.size, dtype=np.complex64)
        val[s != 0] = 1 / (T * s[s != 0])
        return val
    return F


def PT1(T, V=1, dB=False):
    """
    Returns the transfer function F(s) = V / (Ts + 1).
    V may be given in dB.
    """
    if dB:
        V = lin(V)
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        return V / (T * s + 1)
    return F


def PT2(omega, D, V=1, dB=False):
    """
    Returns the transfer function
    F(s) = V / ((s/omega)^2 + 2D/omega * s + 1).
    V may be given in dezibel.
    """
    if dB:
        V = lin(V)
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        return V / ((s / omega) ** 2 + (2 * D / omega) * s + 1)
    return F


def PD1(T, V=1, dB=False):
    """
    Returns the transfer function F(s) = V * (Ts + 1).
    V may be given in dezibel.
    """
    if dB:
        V = lin(V)
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        return T * s + 1
    return F


def GENERAL(a=[1], b=[1]):
    """
    Returns the transfer function
    F(s) = (... a[1] * s + a[0]) / (... b[1] * s + b[0])
    """
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        counter = np.poly1d(a)
        denominator = np.poly1d(b)
        return counter(s) / denominator(s)
    return F

# Composite transfer functions

def PROD(functions):
    """
    Returns the product of the given transfer functions
    as a new transfer function.
    """
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        prod = np.ones(s.shape, dtype=np.complex64)
        for F in functions:
            prod *= F(s)
        return prod
    return F


def SUM(functions):
    """
    Returns the sum of the given transfer functions
    as a new transfer function.
    """
    def F(s):
        s = np.asarray(s, dtype=np.complex64)
        sum = np.zeros(s.shape, dtype=np.complex64)
        for F in functions:
            sum += F(s)
        return sum
    return F


def FEEDBACK_LOOP(F1, F2):
    """
    Returns the transfer function of the feedback loop
    F(s) = F1 / (1 + F1 * F2)
    """
    def F(s):
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
        G = N // 2 + 1

        self.steps = []

        for F in functions:
            # Inverse real fourier transform, only takes lower half of spectrum
            spectrum = F(1j * omega)
            spectrum = spectrum[0:G]
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
