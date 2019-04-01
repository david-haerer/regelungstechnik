#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from abc import ABC, abstractmethod


# Transformations

def unify(val, threshold=0.1):
    prefix = ["",
              "k", "M", "G", "T", "P", "E", "Z", "Y",
              "y", "z", "a", "f", "p", "n", "µ", "m"]

    i = 0

    val = np.asarray(val, dtype=np.float64)

    # Multiple
    if val.max() > 1:
        while val.max() >= threshold * 1e3:
            val *= 1e-3
            i += 1

    # Fraction
    elif val.max() < 1:
        while val.max() < threshold:
            val *= 1e3
            i -= 1

    if val.size == 1:
        val = float(val)

    return val, prefix[i]

def lin(val):
    """
    Returns the linear values of the given dB values.
    """
    return 10 ** (val / 20)

def dB(val):
    """
    Returns the absolute value in dB of the given values.
    """
    return 20 * np.log10(np.abs(val))

def angle(val, min, max):
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

def integrate(y, t):
    delta = t[1] - t[0]
    val = delta * y.cumsum()
    return val

def fft(y):
    l = y.size // N + 1
    val = np.fft.fft(y)
    val[l:] = 0
    return val

def ifft(Y):
    l = Y.size // N + 1
    val = np.fft.irfft(Y[:l])
    return val


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


# Basic functions

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


# Elements

class Element(ABC):
    @abstractmethod
    def H(s):
        """
        The transfer function H(s).
        """
        pass

    @abstractmethod
    def h(t):
        """
        The impulse response h(t).
        """
        pass

    @abstractmethod
    def w(t):
        """
        The step response w(t).
        """
        pass

    @abstractmethod
    def __str__(self):
        pass

class P(Element):
    def __init__(self, V=1, dB=False):
        """
        The basic element P (proportional).
        V may be given in dB.
        """
        super().__init__()

        if dB:
            V = lin(V)

        self.V = V

    def H(s):
        """
        The transfer function H(s) = V.
        """
        s = np.asarray(s, dtype=np.complex64)
        val = self.V * np.ones(shape=s.shape, dtype=np.complex64)
        return val

    def h(t):
        """
        The impulse response h(t) = [V, 0, ...]
        """
        t = np.asarray(t, dtype=np.float64)
        val = np.zeros(shape=t.shape, dtype=np.float64)
        val[t = 0] = self.V
        return val

    def w(t):
        """
        The step response w(t) = V
        """
        t = np.asarray(t, dtype=np.float64)
        val = self.V * np.ones(shape=t.shape, dtype=np.float64)
        val[t < 0] = 0
        return val

    def __str__(self):
        val = "Gain P"
        val += " with " + str_V(self.V)
        return val

class I(Element):
    def __init__(self, T=1):
        """
        The basic element I (integrator).
        """
        super().__init__()

        self.T = T

    def H(s):
        """
        The transfer function H(s) = 1 / Ts.
        """
        s = np.asarray(s, dtype=np.complex64)
        val = 1 / (self.T * s)
        return val

    def h(t):
        """
        The impulse response h(t) = V.
        """
        t = np.asarray(t, dtype=np.float64)
        val = 1 / self.T * np.ones(shape=t.shape, dtype=np.float64)
        val[t < 0] = 0
        return val

    def w(t):
        """
        The step response w(t) = 1 / T * t.
        """
        t = np.asarray(t, dtype=np.float64)
        val = 1 / self.T * t
        val[t < 0] = 0
        return val

    def __str__(self):
        val = "Integrator I"
        val += " with " + str_T(self.T)
        return val

class D(Element):
    def __init__(self, T=1):
        """
        The basic element D (differentiator).
        """
        super().__init__()

        self.T = T

    def H(s):
        """
        The transfer function H(s) = T * s.
        """
        s = np.asarray(s, dtype=np.complex64)
        val = T * s
        return val

    def h(t):
        """
        The impulse response h(t) = 0.
        """
        t = np.asarray(t, dtype=np.float64)
        val = np.zeros(shape=t.shape, dtype=np.float64)
        return val

    def w(t):
        """
        The step response w(t) = [1/T, 0, ...].
        """
        t = np.asarray(t, dtype=np.float64)
        val = np.zeros(shape=t.shape, dtype=np.float64)
        val[t = 0] = 1 / self.T
        return val

    def __str__(self):
        val = "Differentiator D"
        val += " with " + str_T(self.T)
        return val

class PT1(Element):
    def __init__(self, T=1, V=1, dB=False):
        """
        The basic element PT1 (low pass of order 1).
        V may be given in dB.
        """
        super().__init__()

        if dB:
            V = lin(V)

        self.T = T
        self.V = V

    def H(s):
        """
        The transfer function H(s) = V / (Ts + 1).
        """
        s = np.asarray(s, dtype=np.complex64)
        val = self.V / (self.T * s + 1)
        return val

    def h(t):
        """
        The impulse response h(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = self.V / self.T * np.exp(-1 * t / self.T)
        val[t < 0] = 0
        return val

    def w(t):
        """
        The step response w(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = self.V * (1 - np.exp(-1 * t / self.T))
        val[t < 0] = 0
        return val

    def __str__(self):
        val = "Low pass PT1 of order 1"
        val += " with " + str_V(self.V)
        val += " and " + str_T(self.T)
        return val

class PT2(Element):
    def __init__(self, omega=1, D=1, V=1, dB=False):
        """
        The basic element PT2 (low pass of order 2).
        V may be given in dB.
        """
        super().__init__()

        if D < 0:
            print("Error: Damping must not be negative.")
            return None

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
        val = self.V / (s / self.omega ** 2 + 2 * self.D / self.omega * s + 1)
        return val

    def h(t):
        """
        The impulse response h(t) of a PT2 at different dampings.
        """
        t = np.asarray(t, dtype=np.float64)

        V = self.V
        omega = self.omega
        D = self.D

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

    def w(t):
        """
        The step response w(t) of a PT2 at different dampings.
        """
        t = np.asarray(t, dtype=np.float64)

        V = self.V
        omega = self.omega
        D = self.D

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

    def __str__(self):
        val = "Low pass PT2 of order 2"
        val += " with " + str_V(self.V)
        val += " and " + str_omega(self.omega)
        val += " and " + str_V(self.D, dB=False, text="D")
        return val

class PD1(Element):
    def __init__(T=1, V=1, dB=False):
        """
        The basic element PD1 (allowance of order 1).
        V may be given in dB.
        """
        super().__init__()

        if dB:
            V = lin(V)

        self.T = T
        self.V = V

    def H(s):
        """
        The transfer function H(s) = V * (Ts + 1).
        """
        s = np.asarray(s, dtype=np.complex64)
        val = self.V * (self.T * s + 1)
        return val

    def h(t):
        """
        The impulse response h(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = self.V * delta(self.t)
        return val

    def w(t):
        """
        The step response w(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = self.V * (self.T * delta(t) + step(t))
        return val

    def __str__(self):
        val = "Allowance PD1 of order 1"
        val += " with " + str_V(self.V)
        val += " and " + str_T(self.T)
        return val

class PROD(Element):
    def __init__(self, elements):
        """
        The composite element PROD.
        It linkes the given elements in seriell.
        """
        super().__init__()
        self.elements = elements

    def H(s):
        """
        The transfer function H(s) = PROD(elements.H).
        """
        s = np.asarray(s, dtype=np.complex64)
        val = np.ones(s.size, dtype=np.complex64)
        for e in self.elements:
            val *= e.H(s)
        return val

    def h(t):
        """
        The impulse response h(t) = CONV(elements.h).
        """
        t = np.asarray(t, dtype=np.float64)
        val = np.delta(t)
        for e in self.elements:
            val = np.convolve(val, e.h(t), 'same')
        return val

    def w(t):
        """
        The step response w(t) = INT(h(t)).
        """
        t = np.asarray(t, dtype=np.float64)
        val = integrate(self.h(t), t)
        return val

    def __str__(self):
        val = ""
        for e in self.elements:
            val += str(e) + "\n"
        return val

class SUM(Element):
    def __init__(self, elements):
        """
        The composite element SUM.
        It linkes the given elements in parallel.
        """
        super().__init__()
        self.elements = elements

    def H(s):
        """
        The transfer function H(s) = PROD(elements.H).
        """
        s = np.asarray(s, dtype=np.complex64)
        val = np.zeros(s.size, dtype=np.complex64)
        for e in self.elements:
            val += e.H(s)
        return val

    def h(t):
        """
        The impulse response h(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = np.zeros(t)
        for e in self.elements:
            val += e.h(t)
        return val

    def w(t):
        """
        The step response w(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = integrate(self.h(t), t)
        return val

    def __str__(self):
        val = ""
        for e in self.elements:
            val += str(e) + "\n"
        return val


# Diagramms

plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["text.usetex"] = True

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

class Diagramm(ABC):
    """
    Abstract base class for bode, impulse- and step-response diagramm.
    """
    def __init__(self, elements, labels, lang="DE"):
        """
        Takes list of elements and their corresponding labels.
        N is the number of samples.
        The supported languages are DE and EN.
        """
        super().__init__()

        lang = lang.upper()
        if lang != "DE" and lang != "EN":
            print("Error: Supported languages are EN and DE.")
            return None

        self.labels = labels
        self.colors = hue_intervall(len(elements), 1, 0.8)
        self.lang = lang

    @abstractmethod
    def plot(self, pick=None):
        """
        Returns figure of the picked elements.
        """
        pass

    def save(self, pick=None, path="", filename="plot.png"):
        """
        Creates and saves diagramm at path/filename.
        """
        fig = self.plot(pick=pick)
        fig.savefig(path + filename)

    def show(self, pick=None):
        """
        Creates and shows diagramm.
        """
        fig = self.plot(pick=pick)
        fig.show()

class BodeDiagram(Diagramm):
    """
    Bode diagramm of an arbitrary number of transfer functions.
    """
    def __init__(self, elements, labels, start, stop, ticks,
        delta_amp=20, delta_phi=45, N=1024, lang="DE"):
        """
        Takes a list of elements with corresponding labels and
        creates a bode diagramm from 10**start to 10**stop
        with a given list of ticks.
        """
        if super().__init__(elements, labels, lang=lang) is None:
            return None

        self.omega = np.logspace(start=start, stop=stop, num=N)

        ticks = np.arange(ticks[0], ticks[1] + 1)
        self.amp_ticks = delta_amp * ticks
        self.phi_ticks = delta_phi * ticks
        phi_min, phi_max = self.phi_ticks[0], self.phi_ticks[-1]

        self.amps = []
        self.phis = []

        for e in elements:
            spectrum = e.H(1j * self.omega)
            self.amps.append(dB(spectrum))
            self.phis.append(angle(spectrum, phi_min, phi_max))

    def plot(self, pick=None):
        """
        Returns matplotlib figure of the bode diagramm.
        """
        fig, ax_amp = plt.subplots(figsize=(8, 4.5), dpi=240)
        fig.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.95)
        ax_phi = ax_amp.twinx()

        ax_amp.set_xscale("log")
        ax_phi.set_xscale("log")

        amps = np.asarray(self.amps)
        phis = np.asarray(self.phis)
        labels = np.asarray(self.labels)
        colors = np.asarray(self.colors)

        canvas = False

        if pick is not None:
            pick = np.asarray(pick, dtype=int)
            if pick.size == 0:
                canvas = True
            amps = amps[pick]
            phis = phis[pick]
            labels = labels[pick]
            colors = colors[pick]

        for amp, phi, label, color in zip(amps, phis, labels, colors):
            amp_label = r"$|$" + label + r"$|$"
            phi_label = r"$\varphi ($" + label + r"$)$"

            ax_amp.plot(self.omega, amp, label=amp_label,
                       color=color, linewidth=2)
            ax_phi.plot(self.omega, phi, label=phi_label,
                        color=color, linewidth=2, linestyle=":")

        ax_amp.set_xlim(self.omega[0], self.omega[-1])

        ax_amp.set_ylim([self.amp_ticks[0], self.amp_ticks[-1]])
        ax_amp.set_yticks(self.amp_ticks)

        ax_phi.set_ylim([self.phi_ticks[0], self.phi_ticks[-1]])
        ax_phi.set_yticks(self.phi_ticks)

        if self.lang == "DE":
            x_label = r"Kreisfrequenz $\omega \ / \ \frac{1}{s}$"
            amp_label r"Betrag $|H(s)| / dB$"
            phi_label = r"Phase $\varphi(H(s)) / °$"
        else:
            x_label = r"Circular Frequency $\omega \ / \ \frac{1}{s}$"
            amp_label r"Amplitude $|H(s)| / dB$"
            phi_label = r"Phase $\varphi(H(s)) / °$"

        ax_amp.set_xlabel(x_label)
        ax_amp.set_ylabel(amp_label)
        ax_phi.set_ylabel(phi_label)

        if not canvas:
            ax_amp.legend(loc="upper left")
            ax_phi.legend(loc="upper right")

        ax_amp.grid(b=True, which="both", axis="both")

        return fig

class StepResponse(object):
    """
    Step response of an arbitrary number of transfer functions.
    """
    def __init__(self, functions, labels, duration, start=0, lang="DE"):
        """
        Sets functions, labels and colors as attributes.
        Computes sample rate, time and omega from duration and N.
        """
        if super().__init__(elements, labels, lang=lang) is None:
            return None

        self.time = np.linspace(start, duration, 1024)

        self.steps = []

        for e in elements:
            self.steps.append(e.w(self.time))

    def plot(self, pick=None, lim=None):
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

        for w, label, color in zip(steps, labels, colors):
            ax.plot(self.time, w, label=label, color=color, linewidth=2)

        ax.set_xlim(self.time[0], self.time[-1])

        if lim is None:
            if canvas:
                lim = [0, 1]
            else:
                v_min = min([step.min() for step in steps])
                v_max = max([step.max() for step in steps])
                lim = [v_min, v_max]

        ax.set_ylim(lim[0], lim[1])

        if self.lang == "DE":
            self.x_label = r"Zeit $t \ / \ s$"
            self.y_label = r"Sprungantwort $w(t)$"
        else:
            self.x_label = r"Time $t \ / \ s$"
            self.y_label = r"Step Response $w(t)$"

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        if not canvas:
            ax.legend(loc="upper right")

        ax.grid(b=True, which="both", axis="both")

        return fig
