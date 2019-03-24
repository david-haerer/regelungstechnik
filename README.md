# Regelungstechnik

Create diagrams for control theory in python.

**WARNING**

This python module is still under development and may change at any time without notice!

## Contents

1. [Getting Started](#getting-started)
2. [Example](#example)
    1. [Creating Transfer Functions](#example1)
    2. [Creating Bode-Diagrams](#example2)
    3. [Creating Step-Responses](#example3)

## Getting Started <a name="getting-started"></a>

Download `regelungstechnik.py` and `main.py`. The example code in `main.py` should give you a good understanding of how the script works.

## Example <a name="example"></a>

### Creating Transfer Functions <a name="example1"></a>

The example code imports `regelungstechnik as rt` and creates transfer functions.

```python
F1 = rt.PT1(T=2e-3, V=0.2)
F2 = rt.PT2(omega=1000, D=0.2)
F = rt.prod([F1, F2])
```

These functions are grouped in a list and a list of corresponding labels is added. The labels can use LaTeX formatting.

```python
functions = [F, F1, F2]

labels = [
    r"$F = F1 \cdot F2$",
    r"$F_1 = PT_1$",
    r"$F_2 = PT_2$"
]
```

### Creating Bode-Diagrams <a name="example2"></a>

Next up, a Bode-Diagram is created and several plots saved.

```python
bode = rt.BodeDiagram(functions, labels, 1.0, 6.0, ticks=range(-7, 3))
bode.save(pick=[], filename="bode_canvas.png")
bode.save(pick=[0], filename="bode_single.png")
bode.save(filename="bode_all.png")
```

![Bode-Diagram as a canvas](images/bode_canvas.png)

An empty Bode-Diagram can be used as a canvas for hand sketches. It only provides marks on the diagram axis.

![Bode-Diagram of one transfer function](images/bode_single.png)

This Bode-Diagram of the transfer function shows the absolute value in dB and the phase in degrees as the dotted line.

![Bode-Diagram of all transfer functions](images/bode_all.png)

This Bode-Diagram shows all the created transfer functions for comparison.

### Creating Step-Responses <a name="example3"></a>

The last piece of code creates a Step-Response and saves several plots.

```python
step = rt.StepResponse(functions, labels, duration=30e-3)
step.save(pick=[], filename="response_canvas.png", v_max=0.225)
step.save(pick=[0], filename="response_single.png", v_max=0.225)
step.save(filename="response_all.png", v_max=1.6)
```

![Step-Response as a canvas](images/response_canvas.png)

An empty Step-Response can be used as a canvas for hand sketches. It only provides marks on the diagram axis.

![Step-Response of one transfer function](images/response_single.png)

This Step-Response of the transfer function shows how the system reacts to the unit jump.

![Step-Response of all transfer functions](images/response_all.png)

This Step-Response shows all the created transfer functions for comparison.
