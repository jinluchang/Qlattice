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

class AmaVal:

    def __init__(self, val = None, corrections = None):
        # should have the following:
        # self.val: sloppy value
        # self.corrections: list [ (val, description_dict,), ... ]
        # description_dict:
        # { source_specification:
        #   (accuracy_level_relative_to_the_basic_accuracy,
        #    probability_of_having_this_accuracy,),
        # }
        # self.corrections should include the sloppy value in its first element.
        self.val = val
        if corrections is None:
            self.corrections = []
        else:
            self.corrections = corrections

    def __repr__(self):
        return f"AmaVal({self.val},{self.corrections})"

###

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

def ama_apply1(f, x):
    if not isinstance(x, AmaVal):
        return f(x)
    elif isinstance(x, AmaVal):
        val = f(x.val)
        corrections = [ (f(v), d,) for v, d in x.corrections ]
        return AmaVal(val, corrections)
    else:
        assert False

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

def ama_apply2(f, x, y):
    if not isinstance(x, AmaVal) and not isinstance(y, AmaVal):
        return f(x, y)
    elif isinstance(x, AmaVal) and not isinstance(y, AmaVal):
        def f1(x1):
            return f(x1, y)
        return ama_apply1(f1, x)
    elif not isinstance(x, AmaVal) and isinstance(y, AmaVal):
        def f1(y1):
            return f(x, y1)
        return ama_apply1(f1, y)
    elif isinstance(x, AmaVal) and isinstance(y, AmaVal):
        val = f(x.val, y.val)
        corrections = []
        for v_x, d_x in x.corrections:
            for v_y, d_y in y.corrections:
                d = merge_description_dict(d_x, d_y)
                if d is not None:
                    corrections.append((f(v_x, v_y), d,))
        return AmaVal(val, corrections)
    else:
        assert False

def ama_extract(x):
    if not isinstance(x, AmaVal):
        return x
    elif isinstance(x, AmaVal):
        val = x.val
        corrections = x.corrections
        assert isinstance(corrections, list)
        assert corrections
        # keys = [ source_specification, ... ]
        keys = list(corrections[0][1].keys())
        def get_level_prob(key):
            s = set([ d[key] for v, d in corrections ])
            return sorted(list(s))
        dict_level_prob = { k: get_level_prob(k) for k in keys }
        for k, v in dict_level_prob.items():
            assert len(v) >= 2
        # dict_level[key] = list of accuracy levels for this key (source_specification)
        dict_level = { k: [ l for l, prob in v ] for k, v in dict_level_prob.items() }
        # dict_prob[key] = list of probability_of_having_this_accuracy
        dict_prob = { k: [ prob for l, prob in v ] for k, v in dict_level_prob.items() }
        # dict_val[(level, ...)] = v
        dict_val = { tuple([ d[k][0] for k in keys ]): v for v, d in corrections }
        def ama_corr(fixed_levels, remaining_keys):
            if not remaining_keys:
                return dict_val[tuple(fixed_levels)]
            else:
                key = remaining_keys[0]
                rest_keys = remaining_keys[1:]
                levels = dict_level[key]
                probs = dict_prob[key]
                vals = [ ama_corr(fixed_levels + [ l, ], rest_keys) for l in levels ]
                corr = vals[0]
                for i in range(1, len(vals)):
                    corr += (vals[i] - vals[i - 1]) / probs[i]
                return corr
        return ama_corr([], keys)
    else:
        assert False

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
