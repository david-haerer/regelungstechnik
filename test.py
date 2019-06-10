import regelungstechnik as rt
import numpy as np


# Plotting

def plot(elements, labels, name, start, stop, ticks, duration, lim):
    print("Bode Diagram")
    bode = rt.BodeDiagram(elements, labels, start=start, stop=stop, ticks=ticks, lang="DE")
    bode.save(path="tutorium/", filename=name + "_bode")

    print("Step Diagram")
    step = rt.StepResponse(elements, labels, duration=duration, lang="DE")
    step.save(path="tutorium/", filename=name + "_step", lim=lim)


# Proportional

def P():
    Pa = rt.P(V=1)
    Pb = rt.P(V=2)
    Pc = rt.P(V=0.5)
    elements = [Pa, Pb, Pc]
    labels = [
        "$V=1$",
        "$V=2$",
        "$V=0,5$"
    ]
    plot(elements, labels, "P", -2, 2, [-2, 2], 10, [0, 2.5])

def PT1():
    PT1a = rt.PT1(T=1)
    PT1b = rt.PT1(T=5)
    PT1c = rt.PT1(T=0.1)

    elements = [PT1a, PT1b, PT1c]
    labels = [
        "$T = 1s$",
        "$T = 5s$",
        "$T = 0,1s$"
    ]

    plot(elements, labels, "PT1", -1, 3, [-3, 1], 10, [0, 1.2])

def PT2():
    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0)

    elements = [PT2a, PT2b, PT2c, PT2d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0$"
    ]
    plot(elements, labels, "PT2", -1, 3, [-5, 2], 5, [0, 2])

def PT3():
    PT1 = rt.PT1(T=0.1)

    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0)

    PT3a = rt.PROD([PT1, PT2a])
    PT3b = rt.PROD([PT1, PT2b])
    PT3c = rt.PROD([PT1, PT2c])
    PT3d = rt.PROD([PT1, PT2d])

    elements = [PT3a, PT3b, PT3c, PT3d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0$"
    ]
    plot(elements, labels, "PT3", -1, 3, [-6, 2], 5, [0, 2])

def all_P():
    P()
    PT1()
    PT2()
    PT3()


# Differentiator

def D():
    Da = rt.D(T=1)
    Db = rt.D(T=5)
    Dc = rt.D(T=0.5)

    elements = [Da, Db, Dc]
    labels = [
        "$T_D = 1s$",
        "$T_D = 5s$",
        "$T_D = 0,5s$"
    ]

    plot(elements, labels, "D", -2, 2, [-3, 4], 10, [0, 2])

def DT1():
    DT1a = rt.Element([1, 0], [1, 1])
    DT1b = rt.Element([5, 0], [5, 1])
    DT1c = rt.Element([0.1, 0], [0.1, 1])

    elements = [DT1a, DT1b, DT1c]
    labels = [
        "$T = 1s$",
        "$T = 5s$",
        "$T = 0,1s$"
    ]

    plot(elements, labels, "DT1", -3, 3, [-2, 3], 10, [0, 1])

def DT2():
    D = rt.D(T=1)

    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0.05)

    DT2a = rt.PROD([D, PT2a])
    DT2b = rt.PROD([D, PT2b])
    DT2c = rt.PROD([D, PT2c])
    DT2d = rt.PROD([D, PT2d])

    elements = [DT2a, DT2b, DT2c, DT2d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0,05$"
    ]

    plot(elements, labels, "DT2", -1, 3, [-3, 3], 5, [-10, 10])

def DT3():
    D = rt.D(T=1)

    PT1 = rt.PT1(T=0.1)

    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0.05)

    DT3a = rt.PROD([D, PT1, PT2a])
    DT3b = rt.PROD([D, PT1, PT2b])
    DT3c = rt.PROD([D, PT1, PT2c])
    DT3d = rt.PROD([D, PT1, PT2d])

    elements = [DT3a, DT3b, DT3c, DT3d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0,05$"
    ]
    plot(elements, labels, "DT3", -1, 3, [-5, 3], 5, [-10, 10])

def all_D():
    D()
    DT1()
    DT2()
    DT3()


# Integrator

def I():
    Ia = rt.I(T=1)
    Ib = rt.I(T=5)
    Ic = rt.I(T=0.5)

    elements = [Ia, Ib, Ic]
    labels = [
        "$T_i = 1s$",
        "$T_i = 5s$",
        "$T_i = 0,5s$"
    ]

    plot(elements, labels, "I", -1, 4, [-7, 2], 5, [0, 5])

def IT1():
    I = rt.I(T=1)

    PT1a = rt.PT1(T=1)
    PT1b = rt.PT1(T=5)
    PT1c = rt.PT1(T=0.1)

    IT1a = rt.PROD([I, PT1a])
    IT1b = rt.PROD([I, PT1b])
    IT1c = rt.PROD([I, PT1c])

    elements = [IT1a, IT1b, IT1c]
    labels = [
        "$T = 1s$",
        "$T = 5s$",
        "$T = 0,1s$"
    ]

    plot(elements, labels, "IT1", -2, 3, [-7, 2], 10, [0, 10])

def IT2():
    I = rt.I(T=1)

    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0.05)

    IT2a = rt.PROD([I, PT2a])
    IT2b = rt.PROD([I, PT2b])
    IT2c = rt.PROD([I, PT2c])
    IT2d = rt.PROD([I, PT2d])

    elements = [IT2a, IT2b, IT2c, IT2d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0,05$"
    ]

    plot(elements, labels, "IT2", -1, 3, [-7, 1], 2, [0, 2])

def IT3():
    I = rt.I(T=1)

    PT1 = rt.PT1(T=0.5)

    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0.05)

    IT3a = rt.PROD([I, PT1, PT2a])
    IT3b = rt.PROD([I, PT1, PT2b])
    IT3c = rt.PROD([I, PT1, PT2c])
    IT3d = rt.PROD([I, PT1, PT2d])

    elements = [IT3a, IT3b, IT3c, IT3d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0,05$"
    ]
    plot(elements, labels, "IT3", -1, 3, [-9, 1], 2, [0, 2])

def all_I():
    I()
    IT1()
    IT2()
    IT3()


# Allowance

def PD():
    print("PD")

    PDa = rt.Element([1,1],[1e-6,1])
    PDb = rt.Element([5,1],[1e-6,1])
    PDc = rt.Element([0.1,1],[1e-6,1])

    elements = [PDa, PDb, PDc]
    labels = [
        "$T_V = 1s$",
        "$T_V = 5s$",
        "$T_V = 0,1s$"
    ]

    plot(elements, labels, "PD", -3, 3, [-1, 4], 2, [0, 5])

def PDT1():
    print("PDT1")

    PD = rt.PD(T=1)

    PT1a = rt.PT1(T=1)
    PT1b = rt.PT1(T=5)
    PT1c = rt.PT1(T=0.1)

    PDT1a = rt.PROD([PD, PT1a])
    PDT1b = rt.PROD([PD, PT1b])
    PDT1c = rt.PROD([PD, PT1c])

    elements = [PDT1a, PDT1b, PDT1c]
    labels = [
        "$T = 1s$",
        "$T = 5s$",
        "$T = 0,1s$"
    ]

    plot(elements, labels, "PDT1", -2, 3, [-2, 2], 10, [0, 10])

def PDT2():
    print("PDT2")

    PD = rt.PD(T=1)

    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0.05)

    PDT2a = rt.PROD([PD, PT2a])
    PDT2b = rt.PROD([PD, PT2b])
    PDT2c = rt.PROD([PD, PT2c])
    PDT2d = rt.PROD([PD, PT2d])

    elements = [PDT2a, PDT2b, PDT2c, PDT2d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0,05$"
    ]

    plot(elements, labels, "PDT2", -1, 3, [-3, 3], 5, [-7.5, 12.5])

def PDT3():
    print("PDT3")

    PD = rt.PD(T=1)

    PT1 = rt.PT1(T=0.1)

    PT2a = rt.PT2(omega=10, D=1)
    PT2b = rt.PT2(omega=10, D=0.7)
    PT2c = rt.PT2(omega=10, D=0.5)
    PT2d = rt.PT2(omega=10, D=0.05)

    PDT3a = rt.PROD([PD, PT1, PT2a])
    PDT3b = rt.PROD([PD, PT1, PT2b])
    PDT3c = rt.PROD([PD, PT1, PT2c])
    PDT3d = rt.PROD([PD, PT1, PT2d])

    elements = [PDT3a, PDT3b, PDT3c, PDT3d]
    labels = [
        r"$D=1$",
        r"$D=0,7$",
        r"$D=0,5$",
        r"$D=0,05$"
    ]
    plot(elements, labels, "PDT3", -1, 3, [-5, 3], 5, [-6, 8])

def all_PD():
    PD()
    PDT1()
    PDT2()
    PDT3()


# --- SCRIPT ---

#all_P()
#all_D()
#all_I()
all_PD()
