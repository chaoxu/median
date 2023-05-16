# the idea is to implement the data structure for monge matrix

from pandas import *
from bisect import bisect_left
#from types import Any

# the idea is we have functions f_1,...,f_n
# we want to find two positions a and b, so sum min(f_i(a),f_i(b)) is minimized.


def minimum_fix_a(fs, bs, a):
    # given a, we find for each function, which b would be the first 
    # given a picewiselinear unimodal function, find the dagger transform
    A: list[tuple[PiecewiseLinearUnimodal, PiecewiseLinearUnimodal]]
    A = [(dagger_transform(f), f) for f in fs]

    # sort the function by the minimum point value
    A.sort(key=lambda z: z[1].minimum()[0])

    #print("A")
    #print(A)
    #exit(0)

    # we need the initial value and slope

    value = []

    current_value = sum([f.evaluate(0) for f in fs])
    current_slope = sum([f.slope_difference(0) for f in fs])
    p = 0
    for i in range(len(bs)):
        # this should be the first time!
        while A[p][0] < bs[i]:
            # note current slope and value is still at b[i-1]
            # we want to remove the entire of function f!
            current_slope -= A[p][1].slope(bs[i - 1])
            current_value -= A[p][1].evaluate(bs[i - 1])
            current_value += A[p][1].evaluate(a)  # what it suppose to be - actual value of f[p].evaluate(a)
            p += 1
        # do we get to a breakpoint!
        current_value += current_slope * bs[i] - bs[i - 1]
        current_slope += sum([f.slope_difference(bs[i]) for f in fs if f.is_breakpoint(bs[i])])
        value.append(current_value)

    # we need to return the value, which b got the value, and FINALLY, which ones smaller than that b
    minv, minindex = min([(value[i], i) for i in range(len(value))])

    left = []
    right = []
    for i in range(len(A)):
        if A[p][0] < bs[minindex]:
            left.append(A[p][1])
        else:
            right.append(A[p][1])

    return minv, minindex, left, right


# should return the optimum value with elements in as and bs
def optimum(fs, a_s, bs):
    # base case, when a_s is too small. 
    if len(a_s) <= 2:
        result = []
        for i in range(len(a_s)):
            mvi, mbi, _, _ = minimum_fix_a(fs, bs, a_s[i])
            result.append((mvi, a_s[i], mbi))
        return min(result)

    # the idea is to do this recursive approach
    ma = a_s[len(a_s) / 2]
    mv, mb, left, right = minimum_fix_a(fs, bs, ma)

    f_left = SumOfPiecewiseLinearUnimodal(left)
    f_right = SumOfPiecewiseLinearUnimodal(right)

    left_a_s = [x for x in a_s if x in f_left.bs and x<len(a_s) / 2]
    right_a_s = [x for x in a_s if x in f_right.bs and x>len(a_s) / 2]

    left_bs = [x for x in bs if x in f_left.bs and x<=mb]
    right_bs = [x for x in bs if x in f_right.bs and x>=mb]

    return min([(mv, ma, mb), optimum(left, left_a_s, left_bs), optimum(right, right_a_s, right_bs)])

class SumOfPiecewiseLinearUnimodal:
    # list of functions
    # list of breakpoints
    def __init__(self, fs):
        self.fs = fs
        self.compute_breakpoints()
        self.neg_infinity_slope = sum([f.neg_infinity_slope for f in fs])
        self.start_value = self.evaluate(self.bs[0])

    # return a PiecewiseLinearUnimodal function
    def flattern(self):
        # find the slope difference
        delta = []
        for x in self.bs:
            diff = 0.0
            for f in self.has_breakpoint[x]:
                diff+=f.bs_delta[x]
            delta.append(diff)
        return PiecewiseLinearUnimodal(self.start_value, self.bs, delta, self.neg_infinity_slope)

    def compute_breakpoints(self):
        breakpoints = dict([])
        for f in self.fs:
            for x in f.breakpoints:
                if x not in breakpoints:
                    breakpoints[x] = []
                breakpoints[x].append(f)
        self.has_breakpoint = breakpoints
        self.bs = sorted(list(breakpoints.keys()))

    def optimum(self, a_s, bs):
        # base case, when a_s is too small.
        if len(a_s) <= 1:
            result = [(float('inf'), 0, 0)] # make sure something reasonble when empty
            for i in range(len(a_s)):
                mvi, mbi, _, _ = self.naive_minimum_fix_a(bs, a_s[i])
                result.append((mvi, a_s[i], mbi))
            return min(result)

        # the idea is to do this recursive approach
        mid = len(a_s) // 2
        ma = a_s[mid]
        mv, mb, left, right = self.naive_minimum_fix_a(bs, ma)

        # two cases, in the left, we still need to show what to do with right functions

        if right:
            left.append(SumOfPiecewiseLinearUnimodal(right).flattern())
        if left:
            right.append(SumOfPiecewiseLinearUnimodal(left).flattern())

        # in recursion, left need to add an additional linear function
        cases = [(mv, ma, mb)]
        if left:
            f_left = SumOfPiecewiseLinearUnimodal(left)
            left_a_s = [x for x in a_s if x in f_left.bs and x < a_s[mid]]
            left_bs = [x for x in bs if x in f_left.bs and x <= mb]
            cases.append(f_left.optimum(left_a_s, left_bs))
        # what about right?
        if right:
            f_right = SumOfPiecewiseLinearUnimodal(right)
            right_a_s = [x for x in a_s if x in f_right.bs and x > a_s[mid]]
            right_bs = [x for x in bs if x in f_right.bs and x >= mb]
            cases.append(f_right.optimum(right_a_s, right_bs))
        return min(cases)

    def naive_optimum(self, a_s, bs):
        # print("asbs", a_s,bs)
        minimum = float('inf')
        for x in a_s:
            for y in bs:
                if x<=y and self.min_evaluate(x,y) <= minimum:
                    minimum = min(self.min_evaluate(x,y),minimum)
                    ma = x
                    mb = y
        return minimum, ma, mb

    def naive_minimum_fix_a(self, bs, a):
        minv = self.min_evaluate(a, a)
        minindex = a
        for b in bs:
            if b >= a:
                if self.min_evaluate(a, b)<=minv:
                    minv = min(self.min_evaluate(a, b),minv)
                    minindex = b
        left = []
        right = []
        for f in self.fs:
            if f.evaluate(a) <= f.evaluate(minindex):
                left.append(f)
            else:
                right.append(f)
        return minv, minindex, left, right

    def min_evaluate(self, x, y):
        return sum([min(f.evaluate(x),f.evaluate(y)) for f in self.fs])

    def evaluate(self, x):
        return sum([f.evaluate(x) for f in self.fs])
    #def optimum(self):
    #    return optimum(self.fs, self.bs, self.bs)

    def __str__(self):
        s = ""
        for f in self.fs:
            s += f.__str__() + "\n"
        return self.bs.__str__() + "\n" + s

    def print(self):
        print("breakpoints", self.bs)

    def print_table(self):
        print(DataFrame([[self.evaluate(a,b) for b in self.bs] for a in self.bs]))


class PiecewiseLinearUnimodal:
    # need breakpoints breakpoints
    def __init__(self, start_value, bs, delta, neg_infinity_slope=0.0):
        self.breakpoints = bs  # the list of breakpoints
        self.delta = delta  # change in slope
        self.start_value = start_value  # start value is f(b[0]).
        self.neg_infinity_slope = neg_infinity_slope  # slope at negative infinity
        self.bs_delta = dict()
        for i in range(len(bs)):
            self.bs_delta[bs[i]] = delta[i]


    def slope_difference(self, x):
        if x in self.breakpoints:
            return self.delta[self.breakpoints.index(x)]
        return 0.0

    def evaluate(self, x):
        bs = self.breakpoints
        delta = self.delta
        n = len(bs)
        current_value = self.start_value
        current_slope = self.neg_infinity_slope
        for i in range(len(bs)):
            if x < bs[i]:
                current_value += (x - bs[max(i - 1, 0)]) * current_slope
                return current_value
            else:
                current_value += (bs[i] - bs[max(i - 1, 0)]) * current_slope
                current_slope += delta[i]
        if x >= bs[len(bs) - 1]:
            current_value += (x - bs[n - 1]) * current_slope
        return current_value

    def slope(self, x):
        bs = self.breakpoints
        delta = self.delta
        n = len(bs)
        current_slope = self.neg_infinity_slope
        for i in range(len(bs)):
            if x < bs[i]:
                return current_slope
            else:
                current_slope += delta[i]
        return current_slope

    def slopes(self):
        slopes = [self.neg_infinity_slope]
        for d in self.delta:
            slopes.append(slopes[-1] + d)
        return slopes

    def values(self):
        return [self.evaluate(x) for x in self.breakpoints]

    # the position where minimum happens
    def minimum(self):
        minx = self.breakpoints[0]
        minv = self.start_value
        for b in self.breakpoints:
            if self.evaluate(b) < minv:
                minv = self.evaluate(b)
                minx = b
        return minx, minv

    def print(self):
        print("bs", self.breakpoints)
        print("delta", self.delta)
        print("start_value",self.start_value)
        print("neg_infinity_slope",self.neg_infinity_slope)
        print()


# given a list of elements, find the piecewise linear function
def median_to_piecewise_linear(xs):
    xs = sorted(xs)
    n = len(xs)

    # compute the frequency of each element
    values = []
    current_value, current_count = xs[0], 1
    for i in range(1, n):
        if xs[i] == xs[i - 1]:
            current_count += 1
        else:
            values.append((current_value, current_count))
            current_value = xs[i]
            current_count = 1
    values.append((current_value, current_count))

    # breakpoints
    # assuming all values are positive
    # start value = x_i
    start_value = sum(xs) - values[0][0] * n
    neg_infinity_slope = -float(n)

    delta = []
    bs = []
    for v, c in values:
        delta.append(2.0 * c)
        bs.append(v)

    return PiecewiseLinearUnimodal(start_value, bs, delta, neg_infinity_slope)


def slope_to_slope_difference(slope):
    return [slope[i] - slope[i + 1] for i in range(len(slope) - 1)]


def dagger_transform(f):
    # the idea is we start with two lists and converge to the middle
    breakpoints = list(f.breakpoints)
    # we should always make sure we start with a
    fa = f.evaluate(f.breakpoints[0])
    fb = f.evaluate(f.breakpoints[-1])

    if fa < fb:
        # compute the correct a
        a = (f.start_value - fb)/f.neg_infinity_slope + f.breakpoints[0]
        breakpoints = [a] + breakpoints

    breakpoints = [breakpoints[0] - 1.0] + breakpoints
    i = 0
    j = len(breakpoints) - 1
    print(i, j)
    bs = []
    slope = []
    while i != j:
        fa = f.evaluate(breakpoints[i])
        fb = f.evaluate(breakpoints[j])

        print("a,b,f(a),f(b):",i,j,fa,fb)
        if fa <= fb:
            # we need to find the correct a
            a = breakpoints[i - 1] + (fb - f.evaluate(breakpoints[i - 1])) / f.slope(breakpoints[i - 1])
            # compute the correct a
            bs.append(a)
            slope.append(f.slope(a) / f.slope(breakpoints[j - 1]))
            j -= 1
        else:  # a > b
            # we need to find the correct b
            b = breakpoints[j - 1] + (fb - f.evaluate(breakpoints[j - 1])) / f.slope(breakpoints[j - 1])
            # compute the correct b
            bs.append(breakpoints[i])
            slope.append(f.slope(breakpoints[i]) / f.slope(b))
            i += 1
        print("new bs",bs)
        print("new slope",slope)
    #print("BS")
    #print(bs)
    # start_value is what happens at first point, which is the second point. 
    return PiecewiseLinearUnimodal(f.start_value, bs[1:], slope_to_slope_difference(slope), slope[0])


def median_func(list_of_medians):
    return SumOfPiecewiseLinearUnimodal([median_to_piecewise_linear(xs) for xs in list_of_medians])

def naive_median_solver(list_of_medians):
    fs = median_func(list_of_medians)
    return fs.naive_optimum(fs.bs,fs.bs)[0]

def median_solver(list_of_medians):
    fs = median_func(list_of_medians)
    return fs.optimum(fs.bs,fs.bs)[0]