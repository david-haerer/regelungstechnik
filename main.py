import regelungstechnik as rt

# Example transfer function
# F(s) = V / (Ts + 1) * 1 / (s^2/omega^2 + 2D/omega * s + 1)
# Bode Diagramm from 10e1 /s to 10e6 /s

F_1 = rt.make_PT1_func(T=2e-3, V=0.2)
F_2 = rt.make_PT2_func(omega=1000, D=0.2)
F = rt.make_prod_func([F_1, F_2])

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

functions = [F, F_1, F_2]

labels = [
    r"$F = F_1 \cdot F_2$",
    r"$F_1 = PT_1$",
    r"$F_2 = PT_2$"
]

rt.triple_bode(specs, functions, labels, name="example")
