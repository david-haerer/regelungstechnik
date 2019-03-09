import regelungstechnik as rt

# Transfer function F = V * F_1 * F_2:
# V = 0.2
# F_1 = T_1 * s + 1 with T_1 = 10e-6
# F_2 = 1 / (T_2 * s + 1) with T_2 = 2e-6

V = rt.make_V_func(V=0.2)
F_1 = rt.make_DT1_func(T=10e-6)
F_2 = rt.make_PT1_func(T=2e-6)
F = rt.make_prod_func([V, F_1, F_2])

# Bode diagramm: 
# 10**3 < omega < 10**7
# -40dB < abs < 40dB with marks at every 20dB
# -90° < phi < 90° with marks at every 45°

# Labels are given as raw strings for LaTeX formatting.

# 1. Empty bode diagramm as canvas for hand sketches:

bode = rt.EmpytBode(3.0, 7.0, -40, 40, 20, -90, 90, 45)
bode.save(filename="bode_canvas.png")

# 2. Bode diagramm only of product function F:

functions = [F]
labels = [r"$F = V \cdot \frac{T_1 s + 1}{T_2 s + 1}$"]
bode = rt.Bode(functions, labels, 3.0, 7.0, -40, 40, 20, -90, 90, 45)
bode.save(filename="bode_single.png")

# 3. Bode diagramm of all transfer functions:

functions = [F, V, F_1, F_2]
labels = [r"$F = V \cdot F_1 \cdot F_2$", r"$V$", r"$F_1 = DT_1$", r"$F_2 = PT_1$"]
bode = rt.Bode(functions, labels, 3.0, 7.0, -40, 40, 20, -90, 90, 45)
bode.save(filename="bode_multiple.png")
