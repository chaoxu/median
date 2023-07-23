import numpy as np
from itertools import product, compress
from sklearn.tree import DecisionTreeRegressor
import time
from monge import *
import random

def mae(S):
    D = np.concatenate(S)
    return np.abs(D - np.median(D)).sum()


def mask(a, b):
    return [i for i in compress(a, b)]


def brute_force_solution(S):
    repeat = len(S)
    if repeat > 10:
        raise Exception('too big to enumerate')
    mininum = np.inf
    for selector in product([0, 1], repeat=repeat):
        if any(selector) and not all(selector):
            compliment = [not j for j in selector]
            value = mae(mask(S, selector)) + mae(mask(S, compliment))
            if value < mininum:
                mininum = value
    return mininum


def compute_value(X, m=None):
    if m is None:
        m = np.median(X)
    return np.abs(m - X).sum()


def sklearn_solution(S):
    S = sorted(S, key=np.median)
    X = np.repeat(np.arange(len(S)), list(map(len, S)))[:, np.newaxis]
    y = np.hstack(S)
    model = DecisionTreeRegressor(max_depth=1)
    model.fit(X, y)
    return np.sum(np.abs(y - model.predict(X)))


def median_heuristic(S):
    S = sorted(S, key=np.median)
    return min([compute_value(np.hstack(S[:i])) + compute_value(np.hstack(S[i:])) for i in range(1, len(S))])


def piece_wise_linear(E):
    def f(i):
        return np.abs(E - i).sum()


def eval(a, b, fs, selections):
    def g(a, b):
        total = 0.
        for f in fs:
            total += min(f(selections[a]), f(selections[b]))
        return total

    return g


def smawk(rows, cols, lookup):
    # find row maxima
    if not rows: return {}

    stack = []
    for c in cols:
        while len(stack) >= 1 and \
                lookup(rows[len(stack) - 1], stack[-1]) < lookup(rows[len(stack) - 1], c):
            stack.pop()
        if len(stack) != len(rows):
            stack.append(c)

    cols = stack

    result = smawk([rows[i] for i in range(1, len(rows), 2)], cols, lookup)

    c = 0
    for r in range(0, len(rows), 2):
        row = rows[r]
        if r == len(rows) - 1:
            cc = len(cols) - 1  # if r is last row, search through last col
        else:
            cc = c  # otherwise only until pos of max in row r+1
            target = result[rows[r + 1]]
            while cols[cc] != target:
                cc += 1
        result[row] = max([(lookup(row, cols[x]), -x, cols[x]) \
                           for x in range(c, cc + 1)])[2]
        c = cc

    return result


def smawk_base_solution(S):
    elements = np.unique(np.concatenate(S))
    rows = [i for i in range(len(elements))]
    cols = rows[:]

    def g(i, j):
        total = 0.0
        if i > j:
            return -np.inf
        for E in S:
            E = np.array(E)
            total += min(
                np.abs(E - elements[i]).sum(),
                np.abs(E - elements[j]).sum(),
            )
        return -total

    result = smawk(rows, cols, g)
    final_result = min(-g(i, j) for i, j in result.items())
    return final_result


def naive_solutons(S):
    elements = np.unique(np.concatenate(S))
    rows = [i for i in range(len(elements))]
    cols = rows[:]

    def g(i, j):
        total = 0.0
        if i > j:
            return np.inf
        for E in S:
            E = np.array(E)
            total += min(
                np.abs(E - elements[i]).sum(),
                np.abs(E - elements[j]).sum(),
            )
        return total

    final_result = min(g(i, j) for i, j in product(rows, cols))
    return final_result


def benchmark(S, cases):
    for name, func in cases:
        start = time.time()
        result = func(S)
        duration = time.time() - start
        print(f'method: {name}, duration: {duration}, result: {result}')


def main():
    #S = [[747, 446], [17, 749], [997, 312]]
    S=[[0,1,3,4,5],[0,1,2,3,4]]
    #S=[[4, 1, 6, 7], [1, 3, 4, 0], [1, 0, 5, 2]]
    #S = [[6, 6,9,8,1924, 6], [11, 11,3,5,3,19], [14, 0,614,0,6, 6], [19, 19, 19],[19, 33, 19],[100000, 19, 19],[19, 62, 19]]
    for i in range(3):
        S.append(random.sample(range(0, 5), 3))

    for j in range(1):
        S = [[0, 2, 3], [1, 5, 0]]
        for i in range(40):
            S.append(random.sample(range(0, 999999), 41))
        median_solver(S)
        #if smawk_base_solution(S) != median_solver(S):
        #    print(smawk_base_solution(S))
        #    print(median_solver(S))
        #    print(S)
        exit(0)
    benchmark(S, [
        #['brute_force', brute_force_solution],
        ['smawk', smawk_base_solution],
        #['naive', naive_solutons],
        #['naive_median_solver', naive_median_solver],
        ['median_solver', median_solver],
        #['sklearn', sklearn_solution], # heuristics
        ['median', median_heuristic]
    ])
    print(S)


if __name__ == "__main__":
    main()