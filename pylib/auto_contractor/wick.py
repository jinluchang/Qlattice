import copy

class Op:

    def __init__(self, otype : str):
        self.otype = otype

    def is_commute(self) -> bool:
        return True

    def __repr__(self) -> str:
        return f"Op({self.otype})"

    def __add__(self, other):
        return mk_expr(self) + other

    __radd__ = __add__

    def __mul__(self, other):
        return mk_expr(self) * other

    def __rmul__(self, other):
        return mk_expr(other) * self

    def __neg__(self):
        return mk_expr(-1) * mk_expr(self)

    def __sub__(self, other):
        return mk_expr(self) + mk_expr(-1) * other

    def __rsub__(self, other):
        return mk_expr(other) + mk_expr(-1) * self

    def sort(self) -> None:
        pass

    def isospin_symmetric_limit(self) -> None:
        pass

class Qfield(Op):

    def __init__(self, otype : str, flavor : str, position : str, spin : str, color : str):
        Op.__init__(self, otype)
        self.f = flavor
        self.p = position
        self.s = spin
        self.c = color

    def is_commute(self) -> bool:
        return False

    def __repr__(self) -> str:
        if self.s == "auto" and self.c == "auto":
            return f"{self.otype}({self.f!r},{self.p!r})"
        else:
            return f"{self.otype}({self.f!r},{self.p!r},{self.s!r},{self.c!r})"

    def list(self):
        return [self.otype, self.f, self.p, self.s, self.c]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

class Qv(Qfield):

    # act as d / d Hb

    def __init__(self, f, p, s, c):
        # interface function
        Qfield.__init__(self, "Qv", f, p, s, c)

class Qb(Qfield):

    # act as d / d Hv

    def __init__(self, f, p, s, c):
        # interface function
        Qfield.__init__(self, "Qb", f, p, s, c)

class Hv(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "Hv", f, p, s, c)

class Hb(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "Hb", f, p, s, c)

class SHv(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "SHv", f, p, s, c)

class HbS(Qfield):

    def __init__(self, f, p, s, c):
        Qfield.__init__(self, "HbS", f, p, s, c)

class S(Op):

    # propagator

    def __init__(self, flavor : str, p1 : str, p2 : str, s1 : str = "auto", s2 : str = "auto", c1 : str = "auto", c2 : str = "auto"):
        Op.__init__(self, "S")
        self.f = flavor
        self.p1 = p1
        self.p2 = p2
        self.s1 = s1
        self.s2 = s2
        self.c1 = c1
        self.c2 = c2

    def __repr__(self) -> str:
        if self.s1 == "auto" and self.s2 == "auto" and self.c1 == "auto" and self.c2 == "auto":
            return f"{self.otype}({self.f!r},{self.p1!r},{self.p2!r})"
        else:
            return f"{self.otype}({self.f!r},{self.p1!r},{self.p2!r},{self.s1!r},{self.s2!r},{self.c1!r},{self.c2!r})"

    def list(self):
        return [self.otype, self.f, self.p1, self.p2, self.s1, self.s2, self.c1, self.c2]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

    def isospin_symmetric_limit(self) -> None:
        if self.f in ["u", "d",]:
            self.f = "l"

class G(Op):

    # spin matrix

    # tag = "x", "y", "z", "t", "5" for gamma matrices

    def __init__(self, tag : str, s1 : str = "auto", s2 : str = "auto"):
        Op.__init__(self, "G")
        self.tag = tag
        self.s1 = s1
        self.s2 = s2

    def __repr__(self) -> str:
        if self.s1 == "auto" and self.s2 == "auto":
            return f"{self.otype}({self.tag!r})"
        else:
            return f"{self.otype}({self.tag!r},{self.s1!r},{self.s2!r})"

    def list(self):
        return [self.otype, self.tag, self.s1, self.s2]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

def copy_op_index_auto(op : Op):
    if op.otype not in ["S", "G"]:
        return op
    op = copy.copy(op)
    if op.otype in ["S", "G"]:
        op.s1 = "auto"
        op.s2 = "auto"
    if op.otype == "S":
        op.c1 = "auto"
        op.c2 = "auto"
    return op

class Tr(Op):

    # a collection of ops taking the trace

    def __init__(self, ops : list[Op], tag = None):
        Op.__init__(self, "Tr")
        if tag is not None:
            # do not perform check if tag is set
            self.ops = ops
            self.tag = tag
            return
        for op in ops:
            assert op.is_commute()
            assert op.otype in ["S", "G",]
        s = None
        c = None
        for op in ops + ops:
            if op.otype in ["S", "G"]:
                if s is not None:
                    assert s == op.s1
                s = op.s2
            if op.otype == "S":
                if c is not None:
                    assert c == op.c1
                c = op.c2
        if s is not None and c is not None:
            self.tag = "sc"
        elif s is not None:
            self.tag = "s"
        elif c is not None:
            self.tag = "c"
        else:
            self.tag = ""
        for i, op in enumerate(ops):
            ops[i] = copy_op_index_auto(op)
        self.ops = ops

    def __repr__(self) -> str:
        return f"{self.otype}({self.ops!r},{self.tag!r})"

    def list(self):
        return [self.otype, self.tag, self.ops]

    def __eq__(self, other) -> bool:
        return self.list() == other.list()

    def sort(self):
        ops = self.ops
        if len(ops) > 1:
            self.ops = sorted([ ops[i:] + ops[:i] for i in range(len(ops)) ], key = repr)[0]

    def isospin_symmetric_limit(self) -> None:
        for op in self.ops:
            op.isospin_symmetric_limit()

def pick_one(xs, i):
    return xs[i], xs[:i] + xs[i + 1:]

def check_trace_spin_index(ops : list[Op], s : str):
    count1 = 0
    count2 = 0
    i1 = None
    i2 = None
    for i, op in enumerate(ops):
        if op.otype in ["S", "G",]:
            if op.s1 == s:
                i1 = i
                count1 += 1
            if op.s2 == s:
                i2 = i
                count2 += 1
    return count1 == 1 and count2 == 1, i1, i2

def check_trace_color_index(ops : list[Op], c : str):
    count1 = 0
    count2 = 0
    i1 = None
    i2 = None
    for i, op in enumerate(ops):
        if op.otype in ["S",]:
            if op.c1 == c:
                i1 = i
                count1 += 1
            if op.c2 == c:
                i2 = i
                count2 += 1
    return count1 == 1 and count2 == 1, i1, i2

def check_trace_op(ops : list[Op], op : Op):
    if op.otype not in ["S", "G",]:
        return False
    if op.otype in ["S", "G",]:
        if not check_trace_spin_index(ops, op.s1)[0]:
            return False
        if not check_trace_spin_index(ops, op.s2)[0]:
            return False
    if op.otype in ["S",]:
        if not check_trace_color_index(ops, op.c1)[0]:
            return False
        if not check_trace_color_index(ops, op.c2)[0]:
            return False
    return True

def update_trace_sc(op, s, c):
    if op.otype in ["S", "G",]:
        s = op.s2
    if op.otype in ["S",]:
        c = op.c2
    return s, c

def pick_trace_op(ops : list[Op], s, c):
    for i, op in enumerate(ops):
        if not check_trace_op(ops, op):
            continue
        if op.otype in ["S", "G",]:
            if s is not None and s != op.s1:
                continue
        if op.otype in ["S",]:
            if c is not None and c != op.c1:
                continue
        return i, op
    return None

def find_trace(ops : list[Op]):
    # return None or (Tr(tr_ops), remaining_ops,)
    size = len(ops)
    for i, op in enumerate(ops):
        if not check_trace_op(ops, op):
            continue
        masks = [ False for op in ops ]
        tr_ops = []
        s = None
        c = None
        masks[i] = True
        tr_ops.append(op)
        s, c = update_trace_sc(op, s, c)
        while True:
            p_op = pick_trace_op(ops, s, c)
            if p_op is None:
                break
            i2, op2 = pick_trace_op(ops, s, c)
            if masks[i2]:
                return Tr(tr_ops), [op for i, op in enumerate(ops) if not masks[i] ]
            masks[i2] = True
            tr_ops.append(op2)
            s, c = update_trace_sc(op2, s, c)
    return None

def collect_traces(ops : list[Op]) -> list[Op]:
    trs = []
    while True:
        ft = find_trace(ops)
        if ft is None:
            return trs + ops
        tr, ops = ft
        trs.append(tr)

class Term:

    def __init__(self, ops, coef : complex = 1):
        self.coef = complex(coef)
        self.c_ops = []
        self.a_ops = []
        for op in ops:
            if op.is_commute():
                self.c_ops.append(op)
            else:
                self.a_ops.append(op)

    def check_commute(self) -> bool:
        for op in self.c_ops:
            assert op.is_commute()
        for op in self.a_ops:
            assert not op.is_commute()

    def __imul__(self, factor : complex):
        self.coef *= factor

    def __add__(self, other):
        return mk_expr(self) + other

    def __add__(self, other):
        return mk_expr(self) + other

    __radd__ = __add__

    def __mul__(self, other):
        return mk_expr(self) * other

    def __rmul__(self, other):
        return mk_expr(other) * self

    def __neg__(self):
        return mk_expr(-1) * mk_expr(self)

    def __sub__(self, other):
        return mk_expr(self) + mk_expr(-1) * other

    def __rsub__(self, other):
        return mk_expr(other) + mk_expr(-1) * self

    def __repr__(self) -> str:
        return f"Term({self.c_ops + self.a_ops},{self.coef})"

    def sort(self) -> None:
        # only sort commutable factors
        for op in self.c_ops:
            op.sort()
        self.c_ops.sort(key = repr)

    def collect_traces(self) -> None:
        self.c_ops = collect_traces(self.c_ops)

    def isospin_symmetric_limit(self) -> None:
        for op in self.c_ops:
            op.isospin_symmetric_limit()

def combine_two_terms(t1 : Term, t2 : Term):
    if t1.c_ops == t2.c_ops and t1.a_ops == t2.a_ops:
        coef = t1.coef + t2.coef
        if coef == 0.0:
            return Term([], 0.0)
        else:
            return Term(t1.c_ops + t1.a_ops, coef)
    else:
        return None

class Expr:

    def __init__(self, terms):
        self.terms = terms

    def __imul__(self, factor : complex):
        for term in self.terms:
            term *= factor

    def __add__(self, other):
        other = mk_expr(other)
        return Expr(self.terms + other.terms)

    __radd__ = __add__

    def __mul__(self, other):
        other = mk_expr(other)
        terms = []
        for t1 in self.terms:
            for t2 in other.terms:
                coef = t1.coef * t2.coef
                t = Term(t1.c_ops + t2.c_ops + t1.a_ops + t2.a_ops, coef)
                terms.append(t)
        return Expr(terms)

    def __rmul__(self, other):
        return mk_expr(other) * self

    def __neg__(self):
        return mk_expr(-1) * mk_expr(self)

    def __sub__(self, other):
        return mk_expr(self) + mk_expr(-1) * other

    def __rsub__(self, other):
        return mk_expr(other) + mk_expr(-1) * self

    def __repr__(self) -> str:
        return f"Expr({self.terms})"

    def sort(self) -> None:
        for term in self.terms:
            term.sort()
        self.terms.sort(key = repr)

    def combine_terms(self) -> None:
        self.terms = combine_terms_expr(self).terms
        self.terms = drop_zero_terms(self).terms

    def round(self, ndigit : int = 14) -> None:
        # interface function
        sexpr = copy.deepcopy(self)
        for term in sexpr.terms:
            coef = term.coef
            term.coef = complex(round(coef.real, ndigit), round(coef.imag, ndigit))
        sexpr.terms = drop_zero_terms(sexpr).terms
        return sexpr

    def collect_traces(self) -> None:
        for term in self.terms:
            term.collect_traces()

    def isospin_symmetric_limit(self) -> None:
        for term in self.terms:
            term.isospin_symmetric_limit()

    def simplify(self, *, is_isospin_symmetric_limit : bool = False) -> None:
        # interface function
        if is_isospin_symmetric_limit:
            self.isospin_symmetric_limit()
        self.sort()
        self.combine_terms()
        self.collect_traces()
        self.sort()
        self.combine_terms()

def simplified(expr : Expr, *, is_isospin_symmetric_limit : bool = False) -> Expr:
    # interface function
    sexpr = copy.deepcopy(expr)
    sexpr.simplify(is_isospin_symmetric_limit = is_isospin_symmetric_limit)
    return sexpr

def mk_expr(x) -> Expr:
    if isinstance(x, int) or isinstance(x, float) or isinstance(x, complex):
        return Expr([Term([], x),])
    elif isinstance(x, Op):
        return Expr([Term([x,], 1),])
    elif isinstance(x, Term):
        return Expr([x,])
    elif isinstance(x, Expr):
        return x
    else:
        print(x)
        assert False

def combine_terms_expr(expr : Expr) -> Expr:
    if not expr.terms:
        return expr
    terms = []
    term = expr.terms[0]
    for t in expr.terms[1:]:
        ct = combine_two_terms(t, term)
        if ct is None:
            terms.append(term)
            term = t
        else:
            term = ct
    terms.append(term)
    return Expr(terms)

def drop_zero_terms(expr : Expr) -> Expr:
    terms = []
    for t in expr.terms:
        if t.coef != 0.0:
            terms.append(t)
    return Expr(terms)

def op_derivative_exp(op : Op):
    if op.otype == "Qv":
        return Term([SHv(op.f, op.p, op.s, op.c),], -1)
    elif op.otype == "Qb":
        return Term([HbS(op.f, op.p, op.s, op.c),], 1)
    else:
        return None

def op_derivative_op(op : Op, op1 : Op):
    if op.otype == "Qv" and op1.otype == "HbS" and op.f == op1.f:
        return Term([S(op.f, op.p, op1.p, op.s, op1.s, op.c, op1.c),], 1)
    elif op.otype == "Qb" and op1.otype == "SHv" and op.f == op1.f:
        return Term([S(op.f, op1.p, op.p, op1.s, op.s, op1.c, op.c),], 1)
    else:
        return None

def flip_sign(i : int) -> int:
    if i % 2 == 0:
        return 1
    else:
        return -1

def op_derivative_term(op : Op, term : Term) -> Expr:
    coef = term.coef
    c_ops = term.c_ops
    a_ops = term.a_ops
    terms = []
    for i in range(len(a_ops)):
        op1 = a_ops[i]
        dop1 = op_derivative_op(op, op1)
        if dop1 is not None:
            sign = flip_sign(i)
            terms.append(Term(dop1.c_ops + c_ops + a_ops[:i] + dop1.a_ops + a_ops[i+1:], sign * dop1.coef * coef))
    de = op_derivative_exp(op)
    if de is not None:
        sign = flip_sign(len(a_ops))
        terms.append(Term(de.c_ops + c_ops + a_ops + de.a_ops, sign * de.coef * coef))
    return Expr(terms)

def op_push_term(op : Op, term : Term) -> Expr:
    if op.otype == "Qv" or op.otype == "Qb":
        return op_derivative_term(op, term)
    else:
        coef = term.coef
        c_ops = term.c_ops
        a_ops = term.a_ops
        return Expr([Term([op,] + c_ops, a_ops, coef),])

def op_push_expr(op : Op, expr : Expr) -> Expr:
    terms = []
    for term in expr.terms:
        terms += op_push_term(op, term).terms
    return Expr(terms)

def is_hop(op : Op) -> bool:
    if op.otype == "SHv" or op.otype == "HbS":
        return True
    if op.otype == "Hv" or op.otype == "Hb":
        return True
    return False

def has_hops(term : Term) -> bool:
    for op in term.a_ops:
        if is_hop(op):
            return True
    return False

def remove_hops(expr : Expr) -> Expr:
    terms = []
    for term in expr.terms:
        if not has_hops(term):
            terms.append(term)
    return Expr(terms)

def contract_term(term : Term) -> Expr:
    coef = term.coef
    c_ops = term.c_ops
    expr = Expr([Term(c_ops, coef),])
    a_ops = term.a_ops
    for op in reversed(a_ops):
        expr = op_push_expr(op, expr)
    return remove_hops(expr)

def contract_expr(expr: Expr) -> Expr:
    # interface function
    # does not change expr
    all_terms = []
    for term in expr.terms:
        all_terms += contract_term(term).terms
    return Expr(all_terms)

if __name__ == "__main__":
    expr = Expr([
        Term([
            Qb("d", "x1", "s1", "c1"),
            G("5", "s1", "s2"),
            Qv("u", "x1", "s2", "c1"),
            Qb("u", "x2", "s3", "c2"),
            G("5", "s3", "s4"),
            Qv("d", "x2", "s4", "c2"),
            ], 1),
        ])
    print(expr)
    c_expr = contract_expr(expr)
    c_expr.simplify()
    print(c_expr)
    c_expr_check = Expr([Term([Tr([G('5'), S('u','x1','x2'), G('5'), S('d','x2','x1')],'sc')],-1)]) - c_expr
    print(c_expr_check)
    c_expr.simplify(is_isospin_symmetric_limit = True)
    print(c_expr)
    print(Qb("u", "x2", "s23", "c23") * Qv("u", "x1", "s11", "c11") - Qb("d", "x2", "s23", "c23") * Qv("d", "x1", "s11", "c11"))
