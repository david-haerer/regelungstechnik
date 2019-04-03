import regelungstechnik as rt


# Example transfer function is a product of a PT1 and PT2 transfer function
# F(s) = V / (Ts + 1) * 1 / ((s/omega)^2 + 2D/omega * s + 1)

PT1 = rt.PT1(T=2e-3, V=0.2)
PT2 = rt.PT2(omega=1000, D=0.2)
PT3 = rt.PROD([PT1, PT2])


# Make a list of the transfer functions with corresponding labels

elements = [PT3, PT1, PT2]

labels = [
    r"$H = PT_1 \cdot PT_2$",
    r"$PT_1$",
    r"$PT_2$"
]


# Create a Bode-Diagram and save several plots

bode = rt.BodeDiagramm(elements, labels, start=1.0, stop=5.0, ticks=[-7, 2], lang="EN")
bode.save(pick=[], path="images/", filename="bode_canvas.png")
bode.save(pick=[0], path="images/", filename="bode_single.png")
bode.save(path="images/", filename="bode_all.png")


# Create a Step-Response and save several plots

step = rt.StepResponse(elements, labels, duration=30e-3, lang="EN")
step.save(pick=[], path="images/", filename="response_canvas.png",lim=[0, 0.225])
step.save(pick=[0], path="images/", filename="response_single.png",lim=[0, 0.225])
step.save(path="images/", filename="response_all.png", lim=[0, 1.6])
