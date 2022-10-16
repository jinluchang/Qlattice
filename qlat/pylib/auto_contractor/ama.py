#    Qlattice (https://github.com/waterret/qlattice)
#
#    Copyright (C) 2021
#
#    Author: Luchang Jin (ljin.luchang@gmail.com)
#    Author: Masaaki Tomii
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import copy
import qlat as q

class AmaVal:

    def __init__(self, val = None, corrections = None):
        # should have the following:
        # self.val: sloppy value or None
        # self.corrections: list [ (val, description_dict,), ... ]
        # description_dict:
        # { source_specification:
        #   (accuracy_level_relative_to_the_basic_accuracy,
        #    probability_of_having_this_accuracy,),
        # }
        # self.corrections should include the sloppy value and self.val is optional.
        self.val = val
        if corrections is None:
            self.corrections = []
        else:
            self.corrections = corrections

    def __repr__(self):
        return f"AmaVal({self.val},{self.corrections})"

    def __mul__(self, other):
        return ama_msc_mult(self, other)

    def __rmul__(self, other):
        return ama_msc_mult(other, self)

    def __add__(self, other):
        return ama_msc_add(self, other)

    def __radd__(self, other):
        return ama_msc_add(other, self)

###

@q.timer
def mk_ama_val(val, source_specification, val_list, rel_acc_list, prob_list):
    # source_specification need to be unique for each propagator source to ensure proper AMA correction for final result
    # e.g. source_specification = ("point", (12, 2, 3, 4,),)
    # val is the sloppy result
    # val == val_list[0]
    # None value in val_list will be removed automatically
    assert len(val_list) == len(rel_acc_list)
    assert len(val_list) == len(prob_list)
    if val is not val_list[0]:
        assert val == val_list[0]
    corrections = []
    for val_i, rel_acc_i, prob_i in zip(val_list, rel_acc_list, prob_list):
        if val_i is not None:
            corrections.append((val_i, { source_specification: (rel_acc_i, prob_i), },))
    if len(corrections) == 1:
        return val
    return AmaVal(val, corrections)

@q.timer
def ama_apply1_corrections(f, x):
    assert isinstance(x, AmaVal)
    corrections = [ (f(v), d,) for v, d in x.corrections ]
    return corrections

def ama_apply1(f, x):
    if not isinstance(x, AmaVal):
        return f(x)
    else:
        return AmaVal(None, ama_apply1_corrections(f, x))

def ama_counts(x):
    # counts how many times need to compute the val
    if not isinstance(x, AmaVal):
        return 1
    else:
        return len(x.corrections)

def merge_description_dict(d1, d2):
    sd1 = set(d1)
    sd2 = set(d2)
    common_keys = sd1 & sd2
    for key in common_keys:
        if d1[key] != d2[key]:
            return None
    new_keys = sd2 - sd1
    d = d1.copy()
    for key in new_keys:
        d[key] = d2[key]
    return d

@q.timer
def ama_apply2_ama_val(f, x, y):
    assert isinstance(x, AmaVal)
    assert isinstance(y, AmaVal)
    corrections = []
    for v_x, d_x in x.corrections:
        for v_y, d_y in y.corrections:
            d = merge_description_dict(d_x, d_y)
            if d is not None:
                corrections.append((f(v_x, v_y), d,))
    return AmaVal(None, corrections)

@q.timer
def ama_apply2_r_ama_val(f, x, y):
    assert isinstance(y, AmaVal)
    corrections = [ (f(x, v), d,) for v, d in y.corrections ]
    return AmaVal(None, corrections)

@q.timer
def ama_apply2_l_ama_val(f, x, y):
    assert isinstance(x, AmaVal)
    corrections = [ (f(v, y), d,) for v, d in x.corrections ]
    return AmaVal(None, corrections)

def ama_apply2_r(f, x, y):
    if not isinstance(y, AmaVal):
        return f(x, y)
    else:
        return ama_apply2_r_ama_val(f, x, y)

def ama_apply2_l(f, x, y):
    if not isinstance(x, AmaVal):
        return f(x, y)
    else:
        return ama_apply2_l_ama_val(f, x, y)

def ama_apply2(f, x, y):
    if not isinstance(x, AmaVal) and not isinstance(y, AmaVal):
        return f(x, y)
    elif isinstance(x, AmaVal) and not isinstance(y, AmaVal):
        return ama_apply2_l_ama_val(f, x, y)
    elif not isinstance(x, AmaVal) and isinstance(y, AmaVal):
        return ama_apply2_r_ama_val(f, x, y)
    else:
        return ama_apply2_ama_val(f, x, y)

@q.timer
def ama_list(*args):
    l = len(args)
    assert l >= 0
    def f_add(rs, x):
        return rs + [ x, ]
    res = []
    for x in args:
        res = ama_apply2(f_add, res, x)
    return res

def ama_apply(f, *args):
    res = ama_list(*args)
    def f_list(xs):
        return f(*xs)
    return ama_apply1(f_list, res)

@q.timer
def ama_extract_ama_val(x, *, is_sloppy = False):
    corrections = x.corrections
    assert isinstance(corrections, list)
    assert corrections
    # keys = [ source_specification, ... ]
    keys = list(corrections[0][1].keys())
    def get_level_prob(key):
        s = set([ d[key] for v, d in corrections ])
        return sorted(list(s))
    # dict_level_prob[key] = list of accuracy level and probability_of_having_this_accuracy pairs for this key (source_specification)
    dict_level_prob = { k: get_level_prob(k) for k in keys }
    for k, v in dict_level_prob.items():
        assert len(v) >= 2
    # dict_level[key] = list of accuracy levels for this key (source_specification)
    dict_level = { k: [ l for l, prob in v ] for k, v in dict_level_prob.items() }
    # dict_prob[key] = list of probability_of_having_this_accuracy
    dict_prob = { k: [ prob for l, prob in v ] for k, v in dict_level_prob.items() }
    # dict_val[(level, ...)] = v
    dict_val = { tuple([ d[k][0] for k in keys ]): v for v, d in corrections }
    # initial val is sloppy only
    if is_sloppy or (not keys):
        return copy.copy(dict_val[tuple([ dict_level[key][0] for key in keys ])])
    def ama_corr(fixed_levels, remaining_keys):
        if not remaining_keys:
            return dict_val[tuple(fixed_levels)]
        else:
            key = remaining_keys[0]
            rest_keys = remaining_keys[1:]
            levels = dict_level[key]
            probs = dict_prob[key]
            vals = [ ama_corr(fixed_levels + [ l, ], rest_keys) for l in levels ]
            val = copy.copy(vals[0])
            for i in range(1, len(vals)):
                diff = copy.copy(vals[i])
                diff -= vals[i - 1]
                diff *= 1 / probs[i]
                val += diff
            return val
    return ama_corr([], keys)

def ama_extract(x, *, is_sloppy = False):
    if not isinstance(x, AmaVal):
        return x
    elif isinstance(x, AmaVal):
        return ama_extract_ama_val(x, is_sloppy = is_sloppy)
    else:
        assert False

def ama_msc_mult(x, y):
    def f(x, y):
        return x * y
    return ama_apply2(f, x, y)

def ama_msc_add(x, y):
    def f(x, y):
        return x + y
    return ama_apply2(f, x, y)

if __name__ == "__main__":
    v1 = mk_ama_val(1.0, "x", [ 1.0, 1.01, 1.011, ], [ 0, 1, 2, ], [ 1.0, 0.1, 0.02, ])
    v2 = mk_ama_val(2.0, "y", [ 2.0, 2.01, 2.011, ], [ 0, 1, 2, ], [ 1.0, 0.1, 0.02, ])
    v3 = mk_ama_val(3.0, "z", [ 3.0, None, None, ], [ 0, 1, 2, ], [ 1.0, 0.1, 0.02, ])
    def f(x, y):
        return x * y
    v4 = ama_apply2(f, v1, v2)
    v5 = ama_apply2(f, v1, v3)
    v6 = ama_apply2(f, v2, v2)
    print(ama_extract(v1), v1)
    print(ama_extract(v2), v2)
    print(ama_extract(v3), v3)
    print(ama_extract(v4), v4)
    print(ama_extract(v5), v5)
    print(ama_extract(v6), v6)
