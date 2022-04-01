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
