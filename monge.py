#monge backup
# the idea is to implement the data structure for monge matrix

from pandas import *
import numpy as np
from bisect import bisect_left
#from types import Any

# the idea is we have functions f_1,...,f_n
# we want to find two positions a and b, so sum min(f_i(a),f_i(b)) is minimized.


# we try to move to the next value
# advance only work for entire functions
# we try to move to the next value
# advance only work for entire functions
def advance(bs, i, a, f1s, f2s, value, slope):
    # print()
    # print("i",i)
    # print("value",value)
    # print("slope",slope)
    # print("good", f1s)
    # print("failures",f2s)
    # print()
    if i >= len(bs):
        return []
    if i==0:
        slope = sum([f.neg_infinity_slope for f in f1s])
    # f1s and f2s are all lists, it has the property that

    # first update the value to potentially correct one assuming no function is removed
    if i >= 1:
        value += (bs[i] - bs[i - 1]) * slope
        #print("value change",(bs[i] - bs[i - 1]) * slope)
    else:
        value += sum([f.evaluate(bs[0]) for f in f1s])

    #print("before change value", value)
    # then update the slope to potentially correct one
    for f in f1s:
        if bs[i] in f.bs_delta:
            slope += f.bs_delta[bs[i]]

    # # however, we also need to remove certain issues
    failure = -1
    j = 0
    while j<len(f1s):
        f = f1s[j]
        #if bs[i] > a and bs[i]>f.bound: (STILL CAN'T REPLACE. WHY?)
        if bs[i] > a and f.evaluate(bs[i]) > f.evaluate(a):
            # print("ohno",bs[i],a,f.evaluate(bs[i]),f.evaluate(a))
            # this means f(bs[i-1])<=a<f(bs[i])
            # so we overcounted
            value -= f.evaluate(bs[i])
            value += f.evaluate(a)
            slope -= f.slope(bs[i])
            failure = j
            j+=1
        else:
            break

    #print("failuresaaaaaaaaa",failure)
    # remove the failures
    new_f2s = f1s[:failure + 1]
    new_f1s = f1s[failure + 1:]

    reasonable_output = []
    if bs[i] >= a:
        reasonable_output = [(value, i)]
    return advance(bs, i + 1, a, new_f1s, new_f2s, value, slope)+reasonable_output

# def test_advance(fs,bs,a):
#     A = [(max(f.dagger().evaluate(a), f.minimum()[0]), f) for f in fs]
#     #     # sort the function by the minimum point value
#     A.sort(key=lambda x:x[0])
#     for u,v in A:
#         v.bound = u
#     funcs = [v for (u,v) in A]
#     #for f in funcs:
#     #    print(f.breakpoints, dagger_transform(f).evaluate(a), f.minimum()[0])
#     advance(bs, 0, a, funcs, [], 0.0, 0.0)
#     #print()

def minimum_fix_a(fs, bs, a):
    # given a, we find for each function, which b would be the first
    # given a picewiselinear unimodal function, find the dagger transform
    A: list[tuple[float, PiecewiseLinearUnimodal]]
    A = [(max(f.dagger().evaluate(a), f.minimum()[0]), f) for f in fs]
    # sort the function by the minimum point value
    A.sort(key=lambda x:x[0])
    #for i in range(len(fs)):
    #    print("dagger transform evaluated on a", a, A[i][0])
    # we need the initial value and slope
    for u, v in A:
        v.bound = u
    f1s = [y for (x,y) in A]
    #for f in f1s:
#
    #    #print(f.breakpoints, dagger_transform(f).evaluate(a), f.minimum()[0])
    #    #dagger_transform(f).print()
    results = advance(bs, 0, a, f1s, [], 0.0, 0.0)

    minv, minindex = min([(u, -v) for u, v in results])

    return minv, bs[-minindex]

class SumOfPiecewiseLinearUnimodal:
    # list of functions
    # list of breakpoints
    def __init__(self, fs):
        self.fs = fs
        self.cache_evaluate = dict()
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
        #print("running optimum!",a_s,bs)
        # base case, when a_s is too small.
        if len(a_s) <= 1:
            result = [(float('inf'), 0, 0)] # make sure something reasonble when empty
            for i in range(len(a_s)):
                mvi, mbi = minimum_fix_a(self.fs, bs, a_s[i])
                #mvi, mbi = self.naive_minimum_fix_a(bs, a_s[i])

                #assert (mvi == mvi2)
                #assert (mbi == mbi2)
                result.append((mvi, a_s[i], mbi))
            return min(result)

        # the idea is to do this recursive approach
        mid = len(a_s) // 2
        ma = a_s[mid]
        mv, mb = minimum_fix_a(self.fs, bs, ma)
        #mv, mb = self.naive_minimum_fix_a(bs, ma)
        #if (mv != mv2 or mb != mb2):
        #    print("different:",mv,mv2,mb,mb2)
        #    print(self.fs)
        #    print("bs",bs)
        #    print("a",ma)
        left = []
        right = []
        for f in self.fs:
            if f.evaluate(ma) <= f.evaluate(mb):
                left.append(f)
            else:
                right.append(f)

        #print(left, left2)
        #print(right, right2)
        # two cases, in the left, we still need to show what to do with right functions
        new_left = left
        new_right = right
        if right:
            new_left = left+[SumOfPiecewiseLinearUnimodal(right).flattern()]
        if left:
            new_right = right+[SumOfPiecewiseLinearUnimodal(left).flattern()]

        right = new_right
        left = new_left

        # in recursion, left need to add an additional linear function
        cases = [(mv, ma, mb)]
        if left:
            f_left = SumOfPiecewiseLinearUnimodal(left)
            left_a_s = [x for x in a_s if x in f_left.bs and x < a_s[mid]]
            left_bs = [x for x in bs if x in f_left.bs and x <= mb]
            cases.append(f_left.optimum(left_a_s, bs))
        # what about right?
        if right:
            f_right = SumOfPiecewiseLinearUnimodal(right)
            right_a_s = [x for x in a_s if x in f_right.bs and x > a_s[mid]]
            right_bs = [x for x in bs if x in f_right.bs and x >= mb]
            cases.append(f_right.optimum(right_a_s, bs))
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
        #print("naive_fix",a)
        #self.print_table()
        minv = self.min_evaluate(a, a)
        minindex = a
        for b in bs:
            if b >= a:
                if self.min_evaluate(a, b)<=minv:
                    minv = min(self.min_evaluate(a, b),minv)
                    minindex = b
        #print(minv, minindex)
        return minv, minindex

    def min_evaluate(self, x, y):
        return sum([min(f.evaluate(x),f.evaluate(y)) for f in self.fs])

    def evaluate(self, x):
        if x not in self.cache_evaluate:
            self.cache_evaluate[x] = sum([f.evaluate(x) for f in self.fs])
        return self.cache_evaluate[x]
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
        print(DataFrame([[self.min_evaluate(a,b) for b in self.bs] for a in self.bs]))


class PiecewiseLinearUnimodal:
    # need breakpoints breakpoints
    def __init__(self, start_value, bs, delta, neg_infinity_slope=0.0):
        self.breakpoints = bs  # the list of breakpoints
        self.delta = delta  # change in slope
        self.slope_array = []
        self.value = []
        self.start_value = start_value  # start value is f(b[0]).
        self.neg_infinity_slope = neg_infinity_slope  # slope at negative infinity
        self.bs_delta = dict()
        self.bs_index = dict()
        for i in range(len(bs)):
            self.bs_delta[bs[i]] = delta[i]
            self.bs_index[bs[i]] = i
        self.minx = None
        self.minv = None
        self.stored_dagger = None
        self.cache_evaluate = dict()
        self.__build_values()


        # build numpy array for breakpoints
        self.lookup = np.array(self.breakpoints)

        # self.bound = None

    # find the smallest breakpoint
    def prev_index(self, x):
        return np.searchsorted(self.lookup,x,side='left')

    def __build_values(self):
        bs = self.breakpoints
        delta = self.delta
        n = len(bs)
        current_value = self.start_value
        current_slope = self.neg_infinity_slope
        self.minx = None
        self.minv = float('inf')
        self.slope_array.append(current_slope)
        for i in range(len(bs)):
            current_value += (bs[i] - bs[max(i - 1, 0)]) * current_slope
            current_slope += delta[i]
            self.value.append(current_value)
            self.slope_array.append(current_slope)
            if current_value<self.minv:
                self.minx = bs[i]
                self.minv = current_value


    def dagger(self):
        if not self.stored_dagger:
            self.stored_dagger = dagger_transform(self)
        return self.stored_dagger

    def is_breakpoint(self, b):
        return b in self.bs_delta
    def slope_difference(self, x):
        if x in self.breakpoints:
            return self.delta[self.breakpoints.index(x)]
        return 0.0

    def evaluate(self, x):
        if x in self.bs_index:
            return self.value[self.bs_index[x]]

        if x not in self.cache_evaluate:
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
            self.cache_evaluate[x] = current_value
        return self.cache_evaluate[x]

    def slope(self, x):
        if x in self.bs_index:
            return self.slope_array[self.bs_index[x]+1]
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
        #if True:
        if not self.minx:
            minx = self.breakpoints[0]
            minv = self.start_value
            for b in self.breakpoints:
                if self.evaluate(b) < minv:
                    minv = self.evaluate(b)
                    minx = b
            self.minx = minx
            self.minv = minv
        return self.minx, self.minv

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
    return [-(slope[i] - slope[i + 1]) for i in range(len(slope) - 1)]


def dagger_transform(f):
    # the idea is we start with two lists and converge to the middle
    breakpoints = list(f.breakpoints)
    # we should always make sure we start with a
    fa = f.evaluate(f.breakpoints[0])
    fb = f.evaluate(f.breakpoints[-1])
    #print(breakpoints)
    if fa < fb:
        # compute the correct a
        a = f.breakpoints[0]-(f.start_value - fb)/f.neg_infinity_slope
        breakpoints = [a] + breakpoints

    breakpoints = [breakpoints[0] - 1.0] + breakpoints
    i = 0
    j = len(breakpoints) - 1
    bs = []
    slope = []
    while i != j:
        fa = f.evaluate(breakpoints[i])
        fb = f.evaluate(breakpoints[j])

        #print("a,b,f(a),f(b):",i,j,breakpoints[i],breakpoints[j],fa,fb)
        if fa <= fb:
            # we need to find the correct a
            if f.slope(breakpoints[i - 1])!=0:
                a = breakpoints[i - 1] + (fb - f.evaluate(breakpoints[i - 1])) / f.slope(breakpoints[i - 1])
            else:
                a = breakpoints[i - 1]
            if not bs:
                start_value = breakpoints[j]
            # compute the correct a
            bs.append(a)
            if f.slope(breakpoints[j - 1]) != 0:
                slope.append(f.slope(a)/f.slope(breakpoints[j - 1]))
            else:
                slope.append(-1.0)
            j -= 1
        else:  # a > b
            if f.slope(breakpoints[j]) != 0:
                #print("fa,fb",fa,fb)
                b = breakpoints[j] + (fa - fb) / f.slope(breakpoints[j])
                #print("b and sum",b,breakpoints[j],(fa - fb) / f.slope(breakpoints[j]))
            else:
                b = breakpoints[j]
            if not bs:
                start_value = b
            bs.append(breakpoints[i])

            if f.slope(b)!=0:
                slope.append(f.slope(breakpoints[i])/f.slope(b))
            else:
                slope.append(-1.0)
            i += 1
    bs, delta = simplify(bs[1:],slope_to_slope_difference(slope))
    return PiecewiseLinearUnimodal(start_value, bs, delta, slope[0])


# the problem is we have to handle the zero case:

def simplify(bs,delta):
    new_bs=[bs[0]]
    new_delta=[delta[0]]
    for i in range(1,len(bs)):
        if delta[i] != 0:
            new_bs.append(bs[i])
            new_delta.append(delta[i])
    return new_bs, new_delta

def median_func(list_of_medians):
    return SumOfPiecewiseLinearUnimodal([median_to_piecewise_linear(xs) for xs in list_of_medians])

def naive_median_solver(list_of_medians):
    fs = median_func(list_of_medians)
    return fs.naive_optimum(fs.bs,fs.bs)[0]

def median_solver(list_of_medians):
    fs = median_func(list_of_medians)
    return fs.optimum(fs.bs,fs.bs)[0]