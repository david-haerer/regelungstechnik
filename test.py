import regelungstechnik as rt
import numpy as np


def plot(elements, labels, name, start, stop, ticks, duration, lim):
    print("Bode Diagram")
    bode = rt.BodeDiagram(elements, labels, start=start, stop=stop, ticks=ticks, lang="DE")
    bode.save(path="tutorium/", filename=name + "_bode")

    print("Step Diagram")
    step = rt.StepResponse(elements, labels, duration=duration, lang="DE")
    step.save(path="tutorium/", filename=name + "_step", lim=lim)


def P():
    Pa = rt.P(V=1)
    Pb = rt.P(V=2)
    Pc = rt.P(V=0.5)
    elements = [Pa, Pb, Pc]
    labels = ["$V=1$", "$V=2$", "$V=0,5$"]
    plot(elements, labels, "P", -2, 2, [-2, 2], 10, [0, 2.5])

def I():
    Ia = rt.I(T=1)
    Ib = rt.I(T=5)
    Ic = rt.I(T=0.5)
    elements = [Ia, Ib, Ic]
    labels = ["$T_i = 1s$", "$T_i = 5s$", "$T_i = 0,5s$"]
    plot(elements, labels, "I", -1, 4, [-7, 2], 5, [0, 5])

def D():
    Da = rt.D(T=1)
    Db = rt.D(T=5)
    Dc = rt.D(T=0.5)
    elements = [Da, Db, Dc]
    labels = ["$T_D = 1s$", "$T_D = 5s$", "$T_D = 0,5s$"]
    plot(elements, labels, "D", -2, 2, [-5, 4], 10, [0, 2])

def PT1():
    PT1a = rt.PT1(T=1)
    PT1b = rt.PT1(T=5)
    PT1c = rt.PT1(T=0.1)
    elements = [PT1a, PT1b, PT1c]
    labels = ["$T = 1s$", "$T = 5s$", "$T = 0,1s$"]
    plot(elements, labels, "PT1", -1, 3, [-3, 1], 10, [0, 1.5])

def PT2():
    PT1a = rt.PT1(T=0.5)
    PT1b = rt.PT1(T=1)
    PT2a = rt.PROD([PT1a, PT1b])
    PT2b = rt.PT2(omega=10, D=1)
    PT2c = rt.PT2(omega=10, D=0.7)
    PT2d = rt.PT2(omega=10, D=0.5)
    PT2e = rt.PT2(omega=10, D=0.05)
    elements = [PT2a, PT2b, PT2c, PT2d, PT2e]
    labels = [
        r"$T_1 = 0,5s, \ T_2 = 1s$",
        r"$\omega_0 = 10\frac{1}{s}, D=1$",
        r"$\omega_0 = 10\frac{1}{s}, D=0,7$",
        r"$\omega_0 = 10\frac{1}{s}, D=0,5$",
        r"$\omega_0 = 10\frac{1}{s}, D=0,05$"
    ]
    plot(elements, labels, "PT2", -1, 3, [-4, 1], 5, [0, 2])

def PT3():
    PT1a = rt.PT1(T=0.1)
    PT1b = rt.PT1(T=0.5)
    PT1c = rt.PT1(T=1)

    PT2a = rt.PROD([PT1b, PT1c])
    PT2b = rt.PT2(omega=10, D=1)
    PT2c = rt.PT2(omega=10, D=0.7)
    PT2d = rt.PT2(omega=10, D=0.5)
    PT2e = rt.PT2(omega=10, D=0.05)

    PT3a = rt.PROD([PT1a, PT2a])
    PT3b = rt.PROD([PT1a, PT2b])
    PT3c = rt.PROD([PT1a, PT2c])
    PT3d = rt.PROD([PT1a, PT2d])
    PT3e = rt.PROD([PT1a, PT2e])

    elements = [PT3a, PT3b, PT3c, PT3d, PT3e]
    labels = [
        r"$T_1 = 0,1s, \ T_2 = 0,5s, \ T_3 = 1s$",
        r"$T = 0,1s, \ \omega_0 = 10\frac{1}{s}, D=1$",
        r"$T = 0,1s, \ \omega_0 = 10\frac{1}{s}, D=0,7$",
        r"$T = 0,1s, \ \omega_0 = 10\frac{1}{s}, D=0,5$",
        r"$T = 0,1s, \ \omega_0 = 10\frac{1}{s}, D=0,05$"
    ]
    plot(elements, labels, "PT3", -1, 3, [-6, 1], 5, [0, 2])

def DT1():
    DT1a = rt.PROD([D, PT1a])
    DT1b = rt.PROD([D, PT1b])
    DT1c = rt.PROD([D, PT1c])

def DT2():
    DT2a = rt.PROD([D, PT2a])
    DT2b = rt.PROD([D, PT2b])
    DT2c = rt.PROD([D, PT2c])
    DT2d = rt.PROD([D, PT2d])
    DT2e = rt.PROD([D, PT2e])

def DT3():
    DT3a = rt.PROD([D, PT3a])
    DT3b = rt.PROD([D, PT3b])
    DT3c = rt.PROD([D, PT3c])
    DT3d = rt.PROD([D, PT3d])
    DT3e = rt.PROD([D, PT3e])

def PD():
    PD = rt.PD(T=1)

def PDT1():
    PDT1a = rt.PROD([PD, PT1a])
    PDT1b = rt.PROD([PD, PT1b])
    PDT1c = rt.PROD([PD, PT1c])

def PDT2():
    PDT2a = rt.PROD([PD, PT2a])
    PDT2b = rt.PROD([PD, PT2b])
    PDT2c = rt.PROD([PD, PT2c])
    PDT2d = rt.PROD([PD, PT2d])
    PDT2e = rt.PROD([PD, PT2e])

def PDT3():
    PDT3a = rt.PROD([PD, PT3a])
    PDT3b = rt.PROD([PD, PT3b])
    PDT3c = rt.PROD([PD, PT3c])
    PDT3d = rt.PROD([PD, PT3d])
    PDT3e = rt.PROD([PD, PT3e])

def PI():
    PI = rt.PI(T=1)

def PID():
    PIDa = rt.PI(TN=2, TV=1)
    PIDb = rt.PI(TN=2, TV=1, Tv=0.1)


# --- SCRIPT ---

D()
