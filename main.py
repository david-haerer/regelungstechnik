import regelungstechnik as rt


# Example transfer function is a product of a PT1 and PT2 transfer function
# F(s) = V / (Ts + 1) * 1 / ((s/omega)^2 + 2D/omega * s + 1)

F1 = rt.PT1(T=2e-3, V=0.2)
F2 = rt.PT2(omega=1000, D=0.2)
F = rt.prod([F1, F2])


# Make a list of the transfer functions with corresponding labels

functions = [F, F1, F2]

labels = [
    r"$F = F1 \cdot F2$",
    r"$F_1 = PT_1$",
    r"$F_2 = PT_2$"
]


# Create a Bode-Diagram and save several plots

bode = rt.BodeDiagram(functions, labels, start=1.0, stop=6.0, ticks=range(-7, 3))
bode.save(pick=[], path="images/", filename="bode_canvas.png")
bode.save(pick=[0], path="images/", filename="bode_single.png")
bode.save(path="images/", filename="bode_all.png")


# Create a Step-Response and save several plots

step = rt.StepResponse(functions, labels, duration=30e-3)
step.save(pick=[], path="images/", filename="response_canvas.png", v_max=0.225)
step.save(pick=[0], path="images/", filename="response_single.png", v_max=0.225)
step.save(path="images/", filename="response_all.png", v_max=1.6)
