import regelungstechnik as rt

# Example transfer function
# F(s) = V / (Ts + 1) * 1 / (s^2/omega^2 + 2D/omega * s + 1)
# Bode Diagramm from 10e1 /s to 10e6 /s

F_1 = rt.PT1(T=2e-3, V=0.2)
F_2 = rt.PT2(omega=1000, D=0.2)
F = rt.prod([F_1, F_2])

functions = [F, F_1, F_2]

labels = [
    r"$F = F_1 \cdot F_2$",
    r"$F_1 = PT_1$",
    r"$F_2 = PT_2$"
]

bode = rt.BodeDiagram(functions, labels, 1.0, 6.0, ticks=range(-7, 3))
bode.save(pick=[], filename="bode_canvas.png")
bode.save(pick=[0], filename="bode_single.png")
bode.save(filename="bode_all.png")

step = rt.StepResponse(functions, labels, duration=30e-3)
step.save(pick=[], filename="response_canvas.png", v_max=0.225)
step.save(pick=[0], filename="response_single.png", v_max=0.225)
step.save(filename="response_all.png", v_max=1.6)
