from monge import *

### test
# breakpoints [6,11,14,19]



#S = [[1,3,9]]
#, [14, 614, 6], [19, 19, 19],[19, 33, 19],[100000, 19, 19],[19, 62, 19]]


# answer should be 22

S = [[2, 3, 4], [0, 1, 8], [5, 6, 7]]
#median_solver(S)
fs = median_func(S)


fs.print_table()
#print(fs.optimum(fs.bs,fs.bs))
#print(fs.naive_optimum(fs.bs,fs.bs))
#print(minimum_fix_a(fs.fs,fs.bs,1))
#print(fs.naive_minimum_fix_a(fs.bs,1))

f = fs.fs[0]
#for f in fs.fs:
f.print()
df = dagger_transform(f)
df.print()
#    print(df.evaluate(1))

#print(minimum_fix_a(fs.fs, fs.bs, 1))


#print(fs.naive_minimum_fix_a(fs.bs, 1))

#for a in range(6):
#    print()

#mv, mb = minimum_fix_a(fs.fs, fs.bs, 747)
#mv2, mb2 = fs.naive_minimum_fix_a(fs.bs, 747)
#print(mv,mb)
#print(mv2,mb2)

#print(fs.bs)


#df = dagger_transform(f)
#f.print()
#print("")
#df.print()
#exit(0)

#f = median_to_piecewise_linear([1,3,9])
#test_advance([f],f.breakpoints,3)
#f = median_func(S)
#f.print_table()
#print("baseline")


#f.fs[0].print()
#print(f.naive_minimum_fix_a(f.bs,9))

#print("testing")
#print(minimum_fix_a(f.fs,f.bs,9))

#g = f.flattern()
#for x in f.bs:
#    print(f.evaluate(x), g.evaluate(x))
#print(f.naive_minimum_fix_a(f.bs, 11.0))

#f.print_table()
#print(minimum_fix_a(f.fs,f.bs,11.0))
exit(0)
#print(median_solver(S))





F1 = PiecewiseLinearUnimodal(0.0, [0.0, 1.0, 2.0, 3.0, 4.0], [1, 1, 2, 1, 1], -3)
F2 = PiecewiseLinearUnimodal(0.0, [0.0, 1.0, 2.0], [1, 3, 1], -2)

# F1.print()
assert (F1.slopes() == [-3, -2, -1, 1, 2, 3])
assert (F1.values() == [0.0, -2.0, -3.0, -2.0, 0.0])

# F2.print()
assert (F2.slopes() == [-2, -1, 2, 3])
assert (F2.values() == [0.0, -1.0, 1.0])
F1.print()

# test dagger_transform
F3 = dagger_transform(F1)
print(F3.slopes())
print(F3.values())

# next to test SumOfPiecewiseLinearUnimodal
Z = SumOfPiecewiseLinearUnimodal([F1, F1])
Z.print()

# next to test: minimum_fix_a

minimum_fix_a(Z.fs, Z.bs, 3.0)
