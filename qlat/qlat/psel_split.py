import numpy as np

class q:
    from qlat_utils import (
            timer,
            Coordinate,
            smod_coordinate,
            get_fname,
            RngState,
            parallel_map,
            random_permute,
            displayln_info,
            get_num_node,
            )
    from .field_selection import (
            PointsSelection,
            )
    from .mpi_utils import (
            get_mpi_chunk,
            get_comm,
            )

def point_dis_sqr(x, y, total_site):
    return q.smod_coordinate(x - y, total_site).sqr()

class PointsDistanceSet:

    """
    self.total_site
    self.tree
    """

    def __init__(self, psel=None):
        if psel is None:
            self.total_site = None
            self.tree = None
            return
        self.set_psel(psel)

    @q.timer
    def set_psel(self, psel):
        self.total_site = psel.total_site
        self.tree = PointsDistanceTree()
        for xg in psel:
            self.add(xg)

    def add(self, xg):
        self.tree = self.tree.add(xg, self.total_site)

    def check(self):
        if self.total_site is None:
            assert self.tree is None
            return
        assert isinstance(self.total_site, q.Coordinate)
        self.tree.check(self.total_site)

    def count(self):
        return self.tree.count()

    def list(self):
        return self.tree.list()

    def find_closest_point_list(self, xg):
        """
        return mini_dis_sqr, point_list
        #
        IMPORTANT: the point itself, even if it is present in the set, is not included in the search.
        """
        assert isinstance(xg, q.Coordinate)
        return self.tree.find_closest_point_list(xg, self.total_site)

    def find_closest_n_point_list(self, xg, n):
        """
        return point_list
        #
        point_list = [ (dis_sqr, xg,), ... ]
        #
        IMPORTANT: the point itself, even if it is present in the set, is not included in the search.
        """
        assert isinstance(xg, q.Coordinate)
        assert isinstance(n, int)
        return self.tree.find_closest_n_point_list(xg, n, self.total_site)

###

class PointsDistanceTree:

    """
    self.point
    self.range_sqr
    self.tree_list
    """

    def __init__(self):
        self.point = None
        self.range_sqr = None
        self.tree_list = None

    @classmethod
    def mk_with_tree(cls, t0):
        """
        Make a new tree with one more level.
        """
        t = cls()
        assert t0.point is not None
        t.point = t0.point
        if t0.range_sqr is None:
            t.range_sqr = 1
        else:
            t.range_sqr = t0.range_sqr * 4
        t.tree_list = [ t0, ]
        return t

    @classmethod
    def mk_with_point(cls, xg, range_sqr=None):
        """
        Make a leaf for a single point if `range_sqr` is `None` or `0`.
        Otherwise, make a tree containing `xg` with range_sqr at least `range_sqr`.
        """
        t = cls()
        t.point = xg
        t.range_sqr = None
        t.tree_list = None
        if (range_sqr is None) or (range_sqr == 0):
            return t
        t = cls.mk_with_tree(t)
        while t.range_sqr < range_sqr:
            t = cls.mk_with_tree(t)
        return t

    def try_add(self, xg, total_site):
        """
        Try to add the point `xg` to this tree.
        Possibly modify the current tree in place.
        return `True` if `xg` is added in place, return `False` otherwise.
        """
        cls = self.__class__
        if self.point is None:
            return False
        if self.range_sqr is None:
            return False
        dis_sqr = point_dis_sqr(xg, self.point, total_site)
        if dis_sqr > self.range_sqr:
            return False
        if self.tree_list is None:
            self.tree_list = []
        for t in self.tree_list:
            if t.try_add(xg, total_site):
                return True
        t = cls.mk_with_point(xg, self.range_sqr // 4)
        self.tree_list.append(t)
        return True

    def add(self, xg, total_site):
        """
        return new tree with the point `xg` added.
        Possibly modify the current tree in place.
        """
        cls = self.__class__
        if self.point is None:
            return cls.mk_with_point(xg)
        t = self
        while True:
            if t.try_add(xg, total_site):
                return t
            t = cls.mk_with_tree(t)
        assert False

    def check(self, total_site):
        if self.point is None:
            assert self.range_sqr is None
            assert self.tree_list is None
            return
        assert isinstance(self.point, q.Coordinate)
        if self.range_sqr is None:
            assert self.tree_list is None
            return
        assert isinstance(self.range_sqr, int)
        self.range_sqr > 0
        assert isinstance(self.tree_list, list)
        sub_range_sqr = self.range_sqr // 4
        if sub_range_sqr == 0:
            sub_range_sqr = None
        for t in self.tree_list:
            assert t.point is not None
            assert t.range_sqr == sub_range_sqr
            t.check(total_site)
            dis_sqr = point_dis_sqr(self.point, t.point, total_site)
            assert dis_sqr == point_dis_sqr(t.point, self.point, total_site)
            assert dis_sqr <= self.range_sqr

    def count(self):
        if self.point is None:
            return 0
        if self.range_sqr is None:
            return 1
        if self.tree_list is None:
            return 0
        return sum([ t.count() for t in self.tree_list ])

    def list(self):
        l = []
        if self.point is None:
            return l
        if self.range_sqr is None:
            l.append(self.point)
            return l
        l.append(self.range_sqr)
        if self.tree_list is None:
            return l
        assert isinstance(self.tree_list, list)
        if len(self.tree_list) == 0:
            return l
        elif len(self.tree_list) == 1:
            t = self.tree_list[0]
            return t.list()
        for t in self.tree_list:
            v = t.list()
            assert isinstance(v, list)
            while isinstance(v, list) and (len(v) == 1):
                v = v[0]
            l.append(v)
        return l

    def find_closest_point_list(self, xg, total_site):
        """
        return mini_dis_sqr, point_list
        #
        IMPORTANT: the point itself, even if it is present in the set, is not included in the search.
        """
        assert isinstance(xg, q.Coordinate)
        assert isinstance(total_site, q.Coordinate)
        mini_dis_sqr = None
        point_list = []
        pending_list = []
        current_list = [ self, ]
        def consider_point(xg1):
            nonlocal mini_dis_sqr, point_list
            dis_sqr = point_dis_sqr(xg, xg1, total_site)
            if dis_sqr == 0:
                # Excluding the same point
                return
            if (mini_dis_sqr is None) or (mini_dis_sqr > dis_sqr):
                mini_dis_sqr = dis_sqr
                point_list = [ xg1, ]
            elif mini_dis_sqr == dis_sqr:
                point_list.append(xg1)
        def tree_dis_sqr(t):
            if t.point is None:
                return 0
            assert isinstance(t.point, q.Coordinate)
            dis_sqr = point_dis_sqr(xg, t.point, total_site)
            return dis_sqr
        def consider_tree(t):
            nonlocal mini_dis_sqr, point_list, pending_list
            if t.point is None:
                return
            if t.range_sqr is None:
                assert t.tree_list is None
                consider_point(t.point)
                return
            assert isinstance(t.range_sqr, int)
            assert isinstance(t.tree_list, list)
            dis_sqr = point_dis_sqr(xg, t.point, total_site)
            if dis_sqr != 0:
                if (mini_dis_sqr is None) or (mini_dis_sqr > dis_sqr):
                    mini_dis_sqr = dis_sqr
                    point_list = []
                if dis_sqr > mini_dis_sqr + t.range_sqr + 2 * np.sqrt(mini_dis_sqr * t.range_sqr):
                    return
            pending_list += sorted(t.tree_list, key=tree_dis_sqr)
        while len(current_list) > 0:
            for t in current_list:
                consider_tree(t)
            current_list = pending_list
            pending_list = []
        return mini_dis_sqr, point_list

    def find_closest_n_point_list(self, xg, n, total_site):
        """
        return point_list
        #
        point_list = [ (dis_sqr, xg,), ... ]
        #
        IMPORTANT: the point itself, even if it is present in the set, is not included in the search.
        """
        assert isinstance(xg, q.Coordinate)
        assert isinstance(n, int)
        assert isinstance(total_site, q.Coordinate)
        mini_dis_sqr = None
        point_list = []
        def consider_point(xg1):
            nonlocal mini_dis_sqr, point_list
            dis_sqr = point_dis_sqr(xg, xg1, total_site)
            if dis_sqr == 0:
                # Excluding the same point
                return
            if (mini_dis_sqr is None) or (mini_dis_sqr > dis_sqr):
                point_list.append((dis_sqr, xg1,))
                point_list.sort(key=lambda x: x[0])
                if len(point_list) >= n:
                    point_list = point_list[:n]
                    mini_dis_sqr = point_list[-1][0]
                else:
                    assert mini_dis_sqr is None
        def tree_dis_sqr(t):
            if t.point is None:
                return 0
            assert isinstance(t.point, q.Coordinate)
            dis_sqr = point_dis_sqr(xg, t.point, total_site)
            return dis_sqr
        def consider_tree(t):
            if t.point is None:
                return
            if t.range_sqr is None:
                assert t.tree_list is None
                consider_point(t.point)
                return
            assert isinstance(t.range_sqr, int)
            assert isinstance(t.tree_list, list)
            dis_sqr = point_dis_sqr(xg, t.point, total_site)
            if (mini_dis_sqr is not None) and (dis_sqr >= mini_dis_sqr + t.range_sqr + 2 * np.sqrt(mini_dis_sqr * t.range_sqr)):
                return
            pending_list = sorted(t.tree_list, key=tree_dis_sqr)
            for t1 in pending_list:
                consider_tree(t1)
        consider_tree(self)
        if len(point_list) >= n:
            assert mini_dis_sqr == point_list[-1][0]
        else:
            assert mini_dis_sqr is None
        return point_list

###

@q.timer
def find_all_closest_point_list(psel, rs=None, is_parallel=True):
    """
    return all_closest_point_list
    where
    all_closest_point_list = [ (mini_dis_sqr, (xg, [ xg1, ..., ],)), ... ]
    """
    fname = q.get_fname()
    if rs is None:
        rs = q.RngState(f"{fname}")
    assert isinstance(rs, q.RngState)
    pds = PointsDistanceSet(psel)
    def find_closest(xg):
        assert isinstance(xg, q.Coordinate)
        mini_dis_sqr, point_list = pds.find_closest_point_list(xg)
        return (mini_dis_sqr, (xg, point_list,))
    all_closest_point_list = []
    if is_parallel:
        xg_sub_list = q.get_mpi_chunk(list(psel))
        chunksize = max(16, len(psel) // (q.get_num_node() * 128))
        if len(psel) // q.get_num_node() >= 128:
            n_proc = None
        else:
            n_proc = 0
        all_closest_point_sub_list = q.parallel_map(
                find_closest, xg_sub_list,
                chunksize=chunksize, n_proc=n_proc,
                )
        for sub_list in q.get_comm().allgather(all_closest_point_sub_list):
            all_closest_point_list += sub_list
    else:
        for xg in list(psel):
            x = find_closest(xg)
            all_closest_point_list.append(x)
    all_closest_point_list = q.random_permute(all_closest_point_list, rs.split(f"permute"))
    all_closest_point_list.sort(key=lambda x: x[0])
    q.displayln_info(0, f"{fname}: psel.total_site={psel.total_site} ; len(psel)={len(psel)} ; {all_closest_point_list[:5]} .")
    return all_closest_point_list

def find_all_closest_n_point_list_ranking_func_default(dis_sqr_list):
    s = 0.0
    for dis_sqr in dis_sqr_list:
        s += 1 / dis_sqr**3
    if s == 0.0:
        return np.inf
    else:
        return 1 / s**(1/6)

@q.timer(is_timer_fork=True)
def find_all_closest_n_point_list(psel, n, ranking_func=None, rs=None):
    """
    return all_closest_n_point_list
    where
    all_closest_n_point_list = [ (ranking, (xg, [ (dis_sqr1, xg1), ..., ],)), ... ]
    #
    ranking_func(dis_sqr_list) => ranking (the smaller the more closer points)
    """
    fname = q.get_fname()
    if ranking_func is None:
        ranking_func = find_all_closest_n_point_list_ranking_func_default
    if rs is None:
        rs = q.RngState(f"{fname}")
    assert isinstance(rs, q.RngState)
    q.displayln_info(0, f"{fname}: psel.total_site={psel.total_site} ; len(psel)={len(psel)} ; n={n} .")
    pds = PointsDistanceSet(psel)
    def find_closest(xg):
        assert isinstance(xg, q.Coordinate)
        point_list = pds.find_closest_n_point_list(xg, n)
        ranking = ranking_func([ p[0] for p in point_list ])
        return (ranking, (xg, point_list,))
    xg_sub_list = q.get_mpi_chunk(list(psel))
    chunksize = max(16, len(psel) // (q.get_num_node() * 128))
    if len(psel) // q.get_num_node() >= 128:
        n_proc = None
    else:
        n_proc = 0
    all_closest_n_point_sub_list = q.parallel_map(
            find_closest, xg_sub_list,
            chunksize=chunksize, n_proc=n_proc,
            )
    all_closest_n_point_list = []
    for sub_list in q.get_comm().allgather(all_closest_n_point_sub_list):
        all_closest_n_point_list += sub_list
    all_closest_n_point_list = q.random_permute(all_closest_n_point_list, rs.split(f"permute"))
    all_closest_n_point_list.sort(key=lambda x: x[0])
    q.displayln_info(0, f"{fname}: psel.total_site={psel.total_site} ; len(psel)={len(psel)} ; n={n} ; {all_closest_n_point_list[:5]} .")
    return all_closest_n_point_list

@q.timer
def psel_split_that_increase_separation_closest(psel, rs=None):
    """
    split `psel` into `psel1` and `psel2`.
    return psel1, psel2
    """
    assert isinstance(psel, q.PointsSelection)
    fname = q.get_fname()
    if rs is None:
        rs = q.RngState(f"{fname}")
    n_points = len(psel)
    total_site = psel.total_site
    q.displayln_info(0, f"{fname}: total_site={total_site} ; n_points={n_points} .")
    pair_list = find_all_closest_point_list(psel, rs.split(f"find_all_closest_pair_list"))
    xg_set1 = set()
    xg_set2 = set()
    for idx, (mini_dis_sqr, (xg, xg_list,)) in enumerate(pair_list):
        xg = xg.to_tuple()
        xg_list = [ v.to_tuple() for v in xg_list ]
        count1 = 0
        count2 = 0
        for xg1 in xg_list:
            if xg1 in xg_set1:
                count1 += 1
            elif xg1 in xg_set2:
                count2 += 1
        if count1 > count2:
            xg_set2.add(xg)
        elif count2 > count1:
            xg_set1.add(xg)
        else:
            assert count1 == count2
            rsi = rs.split(f"select-{idx}")
            choice = rsi.select([ 1, 2, ])
            if choice == 1:
                xg_set1.add(xg)
            elif choice == 2:
                xg_set2.add(xg)
            else:
                assert False
    xg_list1 = list(xg_set1)
    xg_list2 = list(xg_set2)
    psel1 = q.PointsSelection(total_site, xg_list1)
    psel2 = q.PointsSelection(total_site, xg_list2)
    return psel1, psel2

@q.timer
def psel_split_that_increase_separation_ranking(psel, n, ranking_func=None, rs=None):
    """
    split `psel` into `psel1` and `psel2`.
    return psel1, psel2
    #
    `n` is the number of closest points to be considered for ranking.
    """
    fname = q.get_fname()
    assert isinstance(psel, q.PointsSelection)
    assert isinstance(n, int)
    if ranking_func is None:
        ranking_func = find_all_closest_n_point_list_ranking_func_default
    if rs is None:
        rs = q.RngState(f"{fname}")
    n_points = len(psel)
    total_site = psel.total_site
    q.displayln_info(0, f"{fname}: total_site={total_site} ; n_points={n_points} .")
    all_closest_point_list = find_all_closest_n_point_list(psel, n, ranking_func, rs.split(f"find_all_closest_pair_list"))
    xg_set1 = set()
    xg_set2 = set()
    for idx, (ranking, (xg, point_list,)) in enumerate(all_closest_point_list):
        xg = xg.to_tuple()
        point_list = [ (dis_sqr, xg1.to_tuple(),) for dis_sqr, xg1 in point_list ]
        dis_sqr_list1 = []
        dis_sqr_list2 = []
        for dis_sqr, xg1 in point_list:
            if xg1 in xg_set1:
                dis_sqr_list1.append(dis_sqr)
            elif xg1 in xg_set2:
                dis_sqr_list2.append(dis_sqr)
        ranking1 = ranking_func(dis_sqr_list1)
        ranking2 = ranking_func(dis_sqr_list2)
        if ranking1 > ranking2:
            xg_set1.add(xg)
        elif ranking2 > ranking1:
            xg_set2.add(xg)
        else:
            assert ranking1 == ranking2
            rsi = rs.split(f"select-{idx}")
            choice = rsi.select([ 1, 2, ])
            if choice == 1:
                xg_set1.add(xg)
            elif choice == 2:
                xg_set2.add(xg)
            else:
                assert False
    xg_list1 = list(xg_set1)
    xg_list2 = list(xg_set2)
    psel1 = q.PointsSelection(total_site, xg_list1)
    psel2 = q.PointsSelection(total_site, xg_list2)
    return psel1, psel2

@q.timer(is_timer_fork=True)
def psel_split_that_increase_separation(psel, mode=None, rs=None):
    """
    split `psel` into `psel1` and `psel2`.
    return psel1, psel2
    """
    if mode is None:
        mode = "ranking"
    if mode == "closest":
        return psel_split_that_increase_separation_closest(psel, rs=rs)
    elif mode == "ranking":
        n = 16
        ranking_func = None
        return psel_split_that_increase_separation_ranking(psel, n=n, ranking_func=ranking_func, rs=rs)
    else:
        assert False

@q.timer(is_timer_fork=True)
def find_closest_dis_sqr_for_psel_list(psel_list, is_parallel=True):
    fname = q.get_fname()
    assert isinstance(psel_list, list)
    def find_closest_dis_sqr(psel_idx):
        psel = psel_list[psel_idx]
        assert isinstance(psel, q.PointsSelection)
        if len(psel) == 0:
            return None
        all_closest_point_list = find_all_closest_point_list(psel, is_parallel=not is_parallel)
        assert isinstance(all_closest_point_list, list)
        assert len(all_closest_point_list) > 0
        closest_point_list = all_closest_point_list[0]
        assert isinstance(closest_point_list, tuple)
        assert len(closest_point_list) == 2
        assert isinstance(closest_point_list[0], int)
        assert isinstance(closest_point_list[1], tuple)
        assert len(closest_point_list[1]) == 2
        assert isinstance(closest_point_list[1][0], q.Coordinate)
        assert isinstance(closest_point_list[1][1], list)
        assert len(closest_point_list[1][1]) > 0
        closest_dis_sqr = closest_point_list[0]
        return closest_dis_sqr
    closest_dis_sqr_list = []
    if is_parallel:
        psel_idx_list = list(range(len(psel_list)))
        psel_idx_sub_list = q.get_mpi_chunk(psel_idx_list)
        if len(psel_idx_sub_list) >= 2:
            n_proc = None
        else:
            n_proc = 0
        chunksize = 1
        closest_dis_sqr_sub_list = q.parallel_map(
                find_closest_dis_sqr, psel_idx_sub_list,
                chunksize=chunksize, n_proc=n_proc,
                )
        for sub_list in q.get_comm().allgather(closest_dis_sqr_sub_list):
            closest_dis_sqr_list += sub_list
    else:
        for psel_idx in range(len(psel_list)):
            closest_dis_sqr = find_closest_dis_sqr(psel_idx)
            closest_dis_sqr_list.append(closest_dis_sqr)
    q.displayln_info(0, f"{fname}: {sorted(closest_dis_sqr_list, reverse=True)[:5]}.")
    return closest_dis_sqr_list

@q.timer(is_timer_fork=True)
def psel_split_n_that_increase_separation(psel, num_piece, rs=None):
    """
    return psel_list
    where `len(psel_list) == num_piece`
    """
    assert isinstance(psel, q.PointsSelection)
    assert num_piece >= 1
    fname = q.get_fname()
    if rs is None:
        rs = q.RngState(f"{fname}")
    current_list = [ psel, ]
    pending_list = []
    idx = 0
    while True:
        if len(current_list) + len(pending_list) >= num_piece:
            psel_list = current_list + pending_list
            assert len(psel_list) == num_piece
            find_closest_dis_sqr_for_psel_list(psel_list)
            return psel_list
        if len(current_list) == 0:
            current_list = pending_list
            pending_list = []
        assert len(current_list) > 0
        psel0 = current_list.pop()
        psel1, psel2 = psel_split_that_increase_separation(psel0, mode=None, rs=rs.split(f"{idx}"))
        idx += 1
        pending_list.append(psel1)
        pending_list.append(psel2)
    assert False
