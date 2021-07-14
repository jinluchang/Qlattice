class Op:

    def __init__(self, otype : str):
        self.otype = otype

class Qfield(Op):

    def __init__(self, otype : str, flavor : str, position : str, spin : str, color : str):
        Op.__init__(self, otype)
        self.f = flavor
        self.p = position
        self.s = spin
        self.c = color

class Qv(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "qv", f, p, s, c)

class Qb(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "qb", f, p, s, c)

class Hv(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "yv", f, p, s, c)

class Hb(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "yb", f, p, s, c)

class SHv(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "Syv", f, p, s, c)

class HbS(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "ybS", f, p, s, c)

class Prop(Op):

    def __init__(self, flavor : str, p1 : str, p2 : str, s1 : str, s2 : str, c1 : str, c2 : str):
        Op.__init__(self, "prop")
        self.f = flavor
        self.p1 = p1
        self.p2 = p2
        self.s1 = s1
        self.s2 = s2
        self.c1 = c1
        self.c2 = c2

class Term:

    def __init__(self, coef : complex, ops):
        self.coef = coef
        self.ops = ops

class Expr:

    def __init__(self, terms):
        self.terms = terms

def op_push(op : Op, expr : Expr):



