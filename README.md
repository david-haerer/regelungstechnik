# Regelungstechnik

Create diagramms for control theory in python.

**WARNING**

This python module is still under development and may change at any time without notice!

## How to use it

Download `regelungstechnik.py` and `main.py`. The example code in `main.py` should give you a good understanding of how the script works.

## Example

The example code imports `regelungstechnik as rt` and creates transfer functions. 

```python
F_1 = rt.make_PT1_func(T=2e-3, V=0.2)
F_2 = rt.make_PT2_func(omega=1000, D=0.2)
F = rt.make_prod_func([F_1, F_2])
```

For the plot you have to specify the specs of the bode diagramm in a dictionary.

```python
specs = {
    "start_exp": 1.0,
    "end_exp": 6.0,

    "db_min": -140,
    "db_max": 40,
    "db_delta": 20,

    "phi_min": -315,
    "phi_max": 90,
    "phi_delta": 45
}
```

Then you have to put the transfer functions you want to plot in a list and create a list with corresponding labels. These labels can use LaTeX formatting.

```python
functions = [F, F_1, F_2]

labels = [
    r"$F = F_1 \cdot F_2$",
    r"$F_1 = PT_1$",
    r"$F_2 = PT_2$"
]
```

Then you can create three diffent plots of the transfer function.

```python
rt.triple_bode(specs, functions, labels, name="example")
```

![Canvas bode diagramm](example_canvas.png)

An empty bode diagramm can be used as a canvas for hand sketches. It only provides marks on the diagramm axis.

![Single bode diagramm](example_single.png)

This bode diagramm of the transfer function `F` shows the absolute value in dB and the phase in degrees as the dotted line. 

![Multiple bode diagramm](example_multiple.png)

This bode diagramm additionally shows all the transfer functions `F` is composed of.
