from auto_contractor.wick import *

class Var(Op):

    def __init__(self, name : str):
        Op.__init__(self, "Var")
        self.name = name

    def __repr__(self):
        return f"{self.otype}({self.name!r})"

    def list(self):
        return [self.otype, self.name]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

def add_positions(s, x):
    if isinstance(x, Term):
        for op in x.c_ops:
            add_positions(s, op)
        for op in x.a_ops:
            add_positions(s, op)
    elif isinstance(x, Op):
        if x.otype == "S":
            s.update([x.p1, x.p2])
        elif x.otype == "Tr":
            for op in x.ops:
                add_positions(s, op)
        elif x.otype == "Qfield":
            s.add(x.p)
    elif isinstance(x, Expr):
        for t in x.terms:
            add_positions(s, t)

def add_prop_variables(variables, x):
    prop_counter = 0
    prop_dataset = {} # prop_dataset[prop_repr] = prop_var
    if isinstance(x, list):
        for i, op in enumerate(x):
            if op.otype == "S":
                prop_repr = repr(op)
                if prop_repr in prop_dataset:
                    x[i] = prop_dataset[prop_repr]
                else:
                    while True:
                        prop_counter += 1
                        name = f"S_{prop_counter}"
                        if name not in variables:
                            break
                    variables[name] = op
                    var = Var(name)
                    prop_dataset[prop_repr] = var
                    x[i] = var
            elif op.otype == "Tr":
                add_prop_variables(variables, op.ops)
    elif isinstance(x, Term):
        add_prop_variables(variables, x.c_ops)
    elif isinstance(x, Expr):
        for t in x.terms:
            add_prop_variables(variables, t)

class CExpr:

    def __init__(self, variables, terms, positions = None):
        self.variables = dict(variables)
        self.terms = dict(terms)
        if positions is not None:
            self.positions = positions
        else:
            s = set()
            for t in terms.values():
                add_positions(s, t)
            self.positions = sorted(list(s))

    def __repr__(self) -> str:
        return f"CExpr({sorted(self.variables.items())},{sorted(self.terms.items())},{self.positions})"

    def collect_prop(self):
        variables = self.variables
        for term in self.terms.values():
            add_prop_variables(variables, term)

def mk_cexpr(expr : Expr):
    return CExpr({}, { f"T_{i+1}" : term for i, term in enumerate(expr.terms) })

if __name__ == "__main__":
    expr = Qb("d", "x1", "s1", "c1") * G("5", "s1", "s2") * Qv("u", "x1", "s2", "c1") * Qb("u", "x2", "s3", "c2") * G("5", "s3", "s4") * Qv("d", "x2", "s4", "c2")
    print(expr)
    expr = simplified(contract_expr(expr))
    print(expr.round())
    cexpr = mk_cexpr(expr.round())
    print(cexpr)
    cexpr.collect_prop()
    print(cexpr)
    print(CExpr([('S_1', S('d','x2','x1')), ('S_2', S('u','x1','x2'))],[('T_1', Term([Tr([G('5'), Var('S_1'), G('5'), Var('S_2')],'sc')],(-1+0j)))],['x1', 'x2']))
