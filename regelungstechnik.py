#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from numpy.polynomial import polynomial
from scipy.signal import lti, TransferFunction
import matplotlib.pyplot as plt
import matplotlib.colors as col


# Plotting

plt.rcParams["axes.labelsize"] = 14
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["text.usetex"] = True

def colorspace(num, saturation=1, value=0.8, start=0.0):
    """
    Returns a list of rgb triples visually evenly spaced in regard of hues;
    start markes first color hue in degrees.
    """
    # Matplotlib uses hues in range [0,1].
    start /= 360.0

    hues = np.arange(0.0, 1.0, 1.0 / num)

    # The hues in the hsv color space are visually not evenly distributed.
    # To compensate this effect, we calculate hue**1.5.
    hues = hues ** 1.5

    colors = []
    for hue in hues:
        hsv = ((hue + start) % 1.0, saturation, value)
        colors.append(col.hsv_to_rgb(hsv))

    return colors


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
    val *= step(t)
    return val

def zero(t, dtype=np.float64):
    val = np.zeros(t.size, dtype=dtype)
    return val

def one(t, dtype=np.float64):
    val = np.ones(t.size, dtype=dtype)
    return val


# Transformations

def unify(val, threshold=0.1):
    prefix = [r"",
              r"k", r"M", r"G", r"T", r"P", r"E", r"Z", r"Y",
              r"y", r"z", r"a", r"f", r"p", r"n", r"µ", r"m"]

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

def angle_norm(phi, min, max, convert=False):
    """
    Returns the phase value in degrees in the range between
    v_min and v_max of the given values.
    """
    if convert:
        phi = np.angle(phi, deg=True)
    # Move phi down
    while (phi > max).any():
        phi[phi > max] -= 360
    # Move phi up
    while (phi < min).any():
        phi[phi < min] += 360
    # Move phi down as far as possible
    while (phi - 360 >= min).any():
        phi[phi - 360 > min] -= 360
    return phi

def fft(y):
    l = y.size // 2 + 1
    val = np.fft.fft(y)
    val[l:] = 0
    return val

def ifft(Y):
    l = Y.size // 2 + 1
    val = np.fft.irfft(Y[:l])
    return val

def polyadd(a, b):
    """
    Add polynomial a by polynomial b.
    a and b are lists from highest order term to lowest.
    """
    val = polynomial.polyadd(a[::-1], b[::-1])
    val = val[::-1]
    return val

def polymul(a, b):
    """
    Multiply polynomial a by polynomial b.
    a and b are lists from highest order term to lowest.
    """
    val = polynomial.polymul(a[::-1], b[::-1])
    val = val[::-1]
    return val


# Elements

class Element(object):
    def __init__(self, counter, denominator):
        """
        An element is described by the counter and denominator of its transfer function. It has corresponding scipy lti and TransferFunction object as attributes.
        counter and denominator are list from highest order term to lowest.
        """
        self.counter = np.asarray(counter, dtype=np.float64)
        self.denominator = np.asarray(denominator, dtype=np.float64)
        self.sys = lti(counter, denominator)
        self.tf = TransferFunction(counter, denominator)

    def H(self, s):
        """
        The transfer function H(s).
        """
        s = np.asarray(s, dtype=np.complex64)
        val = self.sys.freqresp()[1]
        return val

    def bode(self, omega):
        """
        The Bode amplitude and phase data.
        """
        omega = np.asarray(omega, dtype=np.complex64)
        omega, amp, phi = self.sys.bode(w=omega, n=512)
        return omega, amp, phi

    def h(self, t):
        """
        The impulse response h(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = self.sys.impulse(T=t, N=512)[1]
        return val

    def w(self, t):
        """
        The step response w(t).
        """
        t = np.asarray(t, dtype=np.float64)
        val = self.sys.step(T=t, N=512)[1]
        return val

    def __add__(self, other):
        """
        self + other
        """
        return SUM([self, other])

    def __radd__(self, other):
        """
        other + self
        """
        return self.__add__(other)

    def __iadd__(self, other):
        """
        self += other
        """
        return self.__add__(other)

    def __mul__(self, other):
        return PROD([self, other])

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imull(self, other):
        return self.__mul__(other)

    def __matmul__(self, other):
        return FEEDBACK(self, other)

    def __rmatmul__(self, other):
        return FEEDBACK(other, self)

    def __imatmul__(self, other):
        return self.__matmul__(other)

class P(Element):
    def __init__(self, V=1, dB=False):
        """
        The basic element P (proportional).
        V may be given in dB.
        """
        if dB:
            V = lin(V)

        self.V = V

        super().__init__([V], [1])

class I(Element):
    def __init__(self, T=1):
        """
        The basic element I (integrator).
        """
        self.T = T

        super().__init__([1], [T, 0])

class D(Element):
    def __init__(self, T=1, real=True):
        """
        The basic element D (differentiator).
        """
        self.T = T

        if real:
            Tv = T * 1e-6
        else:
            Tv = 0

        super().__init__([T, 0], [Tv, 1])

class PT1(Element):
    def __init__(self, T=1, V=1, dB=False):
        """
        The basic element PT1 (low pass of order 1).
        V may be given in dB.
        """
        if dB:
            V = lin(V)

        self.T = T
        self.V = V

        super().__init__([V], [T, 1])

class PT2(Element):
    def __init__(self, omega=1, D=1, V=1, dB=False):
        """
        The basic element PT2 (low pass of order 2).
        V may be given in dB.
        """
        if D < 0:
            print("Error: Damping must not be negative.")
            return None

        if dB:
            V = lin(V)

        self.omega = omega
        self.D = D
        self.V = V

        super().__init__([V], [1 / (omega ** 2), 2 * D / omega, 1])

class PD(Element):
    def __init__(self, T=1, V=1, dB=False, Tv=0):
        """
        The basic element PD1 (allowance of order 1).
        V may be given in dB.
        """
        if dB:
            V = lin(V)

        self.T = T
        self.V = V
        self.Tv = Tv
        
        if Tv == 0:
            self.Tv = T * 1e-6

        super().__init__([T * V, V], [Tv, 1])

class PI(Element):
    def __init__(self, V=1, T=1, db=False):
        if dB:
            V = lin(V)

        self.T = T
        self.V = V

        super().__init__([V * T, V], [T, 0])

class PID(Element):
    def __init__(self, V=1, TN=1, TV=1, db=False, Tv=0):
        if dB:
            V = lin(V)

        self.V = V
        self.TN = TN
        self.TV = TV
        self.Tv = Tv

        super().__init__([V * TN * TV, V * (TN + TV), V], [TN * Tv, TN, 0])

class PROD(Element):
    def __init__(self, elements):
        """
        The composite element PROD.
        It linkes the given elements in seriell.
        """
        self.elements = elements[:]

        # Cancel polynomials

        cs = [e.counter for e in self.elements]
        ds = [e.denominator for e in self.elements]

        for i in range(len(elements)):
            for j in range(len(elements)):
                if cs[i].size == ds[j].size and (cs[i] == ds[j]).all():
                    cs[i] = np.ones(1)
                    ds[j] = np.ones(1)
                    break

        counter, denominator = [1], [1]

        for c in cs:
            counter = polymul(counter, c)

        for d in ds:
            denominator = polymul(denominator, d)

        super().__init__(counter, denominator)

class SUM(Element):
    def __init__(self, elements):
        """
        The composite element SUM.
        It linkes the given elements in parallel.
        """
        self.elements = elements

        counter, denominator = [0], [1]

        for e in elements:
            expander = e.counter
            remains = elements[:].remove(e)
            for r in remains:
                expander = polymul(expander, r.denominator)
            counter = polyadd(counter, expander)
            denominator = polymul(denominator, e.denominator)

        super().__init__(counter, denominator)

class FEEDBACK(Element):
    def __init__(self, feed, back):
        """
        The composite element FEEDBACK = feed / (1 + feed * back)
        """
        counter = polymul(feed.counter, back.denominator)
        denominator_1 = polymul(feed.denominator, back.denominator)
        denominator_2 = polymul(feed.counter, back.counter)
        denominator = polyadd(denominator_1, denominator_2)

        super().__init__(counter, denominator)


# Diagrams

class BodeDiagram(object):
    """
    Bode diagram of an arbitrary number of transfer functions.
    """
    def __init__(self, elements, labels, start, stop, ticks, delta_amp=20, delta_phi=45, N=1024, lang="DE"):
        """
        Takes a list of elements with corresponding labels and
        creates a bode diagram from 10**start to 10**stop
        with a given list of ticks.
        """
        lang = lang.upper()
        if lang != "DE" and lang != "EN":
            print("Error: Supported languages are EN and DE.")
            print("Language is set to EN.")
            lang = "EN"

        self.labels = labels
        self.colors = colorspace(len(elements))
        self.lang = lang

        self.omega = np.logspace(start=start, stop=stop, num=N)

        ticks = np.arange(ticks[0], ticks[1] + 1)
        self.amp_ticks = delta_amp * ticks
        self.phi_ticks = delta_phi * ticks

        self.amps = []
        self.phis = []

        for e in elements:
            omega, amp, phi = e.bode(self.omega)
            self.amps.append(amp)
            self.phis.append(angle_norm(phi, self.phi_ticks[0], self.phi_ticks[-1]))

    def plot(self, pick=None):
        """
        Returns matplotlib figure of the bode diagram.
        """
        fig, ax_amp = plt.subplots(figsize=(8, 4.5), dpi=240)
        fig.subplots_adjust(left=0.125, right=0.875, bottom=0.15, top=0.95)
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
            if pick.size == 0 or len(amps) == 0:
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

        x_label = r"Circular Frequency "
        amp_label = r"Amplitude "
        phi_label = r"Phase "

        if self.lang == "DE":
            x_label = r"Kreisfrequenz "
            amp_label = r"Betrag "
            phi_label = r"Phase "

        x_label += r"$\ \omega \ / \ \frac{1}{s}$"
        amp_label += r"$\ |H(j \omega)| \ / \ dB$"
        phi_label += r"$\ \varphi(H(j \omega)) \ / \ ^\circ$"

        ax_amp.set_xlabel(x_label)
        ax_amp.set_ylabel(amp_label)
        ax_phi.set_ylabel(phi_label)

        if not canvas:
            ax_amp.legend(loc="upper left")
            ax_phi.legend(loc="upper right")

        ax_amp.grid(b=True, which="both", axis="both")

        return fig

    def save(self, pick=None, path="", filename="plot.png"):
        """
        Creates and saves diagram at path/filename.
        """
        fig = self.plot(pick=pick)
        fig.savefig(path + filename)

    def show(self, pick=None):
        """
        Creates and shows diagram.
        """
        fig = self.plot(pick=pick)
        fig.show()

class StepResponse(object):
    """
    Step response of an arbitrary number of transfer functions.
    """
    def __init__(self, elements, labels, duration, start=0, lang="DE"):
        """
        Sets functions, labels and colors as attributes.
        Computes sample rate, time and omega from duration and N.
        """
        lang = lang.upper()
        if lang != "DE" and lang != "EN":
            print("Error: Supported languages are EN and DE.")
            print("Default: Language is set to EN.")
            lang = "EN"

        self.labels = labels
        self.colors = colorspace(len(elements))
        self.lang = lang

        self.time = np.linspace(start, duration, 1024)

        self.steps = []

        for e in elements:
            self.steps.append(e.w(self.time))

        self.time, self.time_prefix = unify(self.time)

    def plot(self, pick=None, lim=None):
        """
        Returns matplotlib figure of the step responses.
        """
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=240)
        fig.subplots_adjust(left=0.125, right=0.925, bottom=0.15, top=0.95)

        steps = np.asarray(self.steps)
        labels = np.asarray(self.labels)
        colors = np.asarray(self.colors)

        canvas = False

        if pick is not None:
            pick = np.asarray(pick, dtype=int)
            if pick.size == 0 or len(steps) == 0:
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

        x_label = r"Time "
        y_label = r"Step Response "

        if self.lang == "DE":
            x_label = r"Zeit "
            y_label = r"Sprungantwort "

        x_label += r"$\ t \ / \ " + self.time_prefix + r"s$"
        y_label += r"$\ w(t)$"

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)

        if not canvas:
            ax.legend(loc="upper right")

        ax.grid(b=True, which="both", axis="both")

        return fig

    def save(self, pick=None, lim=None, path="", filename="plot.png"):
        """
        Creates and saves diagram at path/filename.
        """
        fig = self.plot(pick=pick, lim=lim)
        fig.savefig(path + filename)

    def show(self, pick=None, lim=None):
        """
        Creates and shows diagram.
        """
        fig = self.plot(pick=pick, lim=lim)
        fig.show()


# Digital Controller

class DigitalPID(object):
    def __init__(self, V, Tn, Tv, delta_t, umax=None):
        self.delta_t = delta_t

        self.KP = V * (Tn + Tv) / Tn
        self.KI = V * delta_t / Tn
        self.KD = V * Tv / delta_t

        self.umax = umax
        self.e_old = 0
        self.ui = 0

    def __call__(self, target, real, verbose=False):
        e = target - real

        up = self.KP * e

        ui = self.ui + self.KI * e

        if self.umax is not None:
            ui = min(ui, self.umax)
            ui = max(ui, -self.umax)

        ud = self.KD * (e - self.e_old)

        u = up + ui + ud

        if self.umax is not None:
            u = min(u, self.umax)
            u = max(u, -self.umax)

        self.e_old = e
        self.ui = ui

        if verbose:
            print(f"e={e:.2f} up={up:.2f} ui={ui:.2f} ud={ud:.2f} u={u:.2f}")

        return u

    def w(self, t):
        t = np.asarray(t, dtype=np.float64)
        s = step(t)
        n = t.size
        val = np.zeros(n, dtype=np.float64)

        tick = t[0]
        val[0] = self(s[0], 0)

        for i in range(1, n):
            if t[i] - tick >= self.delta_t:
                tick = t[i]
                val[i] = self(s[i], 0)
            else:
                val[i] = val[i - 1]

        return val

    def code(self, lang="DE"):
        target = "soll" if lang == "DE" else "target"
        real = "ist" if lang == "DE" else "real"
        e_old = "e_alt" if lang == "DE" else "e_old"

        code = f"float KP = {self.KP:.2f};\n"
        code += f"float KI = {self.KI:.2f};\n"
        code += f"float KD = {self.KD:.2f};\n"

        if self.umax is not None:
            code += f"float umax = {self.umax};\n"

        code += f"\n"
        code += f"float regler(float {target}, float {real}) {{\n"
        code += f"    float e, u, up, ud;\n"
        code += f"    static float ui = 0;\n"
        code += f"    static float {e_old} = 0;\n"
        code += f"\n"
        code += f"    up = KP * e;\n"
        code += f"    ui = ui + KI * e;\n"
        code += f"    ud = KD * (e - {e_old});\n"
        code += f"\n"

        if self.umax is not None:
            code += f"    if(ui > umax);\n"
            code += f"        ui = umax;\n"
            code += f"    if(ui < -umax);\n"
            code += f"        ui = -umax;\n"
            code += f"\n"

        code += f"    u = up + ui + ud;\n"
        code += f"\n"

        if self.umax is not None:
            code += f"    if(u > umax);\n"
            code += f"        u = umax;\n"
            code += f"    if(u < -umax);\n"
            code += f"        u = -umax;\n"
            code += f"\n"

        code += f"    {e_old} = e;\n"
        code += f"\n"
        code += f"    return u;\n"
        code += f"}}"

        return code

    def reset(self):
        self.e_old = 0
        self.ui = 0
