from monge import *

### test
# breakpoints [6,11,14,19]



#S = [[1,3,9]]
#, [14, 614, 6], [19, 19, 19],[19, 33, 19],[100000, 19, 19],[19, 62, 19]]


# answer should be 22

S = [[70974, 79595, 15251, 40792, 28026, 54394, 95874, 10114, 66999, 19974, 78746], [94332, 29205, 60427, 38707, 22017, 12287, 29308, 24897, 52921, 6168, 10657], [47661, 66183, 78421, 52249, 62497, 85434, 93775, 31453, 47601, 72508, 15315], [82642, 11733, 83785, 16950, 74766, 87251, 86583, 26090, 45804, 64099, 63586], [48126, 59025, 14587, 61640, 86763, 74582, 66265, 83673, 34623, 63346, 19142], [87436, 12011, 85101, 71730, 59497, 17389, 59558, 86967, 63157, 78130, 87342], [12596, 9861, 53689, 33345, 22761, 81052, 23941, 62362, 78964, 10361, 46864], [89456, 37188, 39876, 45293, 85115, 26745, 98333, 38679, 61766, 37898, 54778], [42503, 75494, 50053, 46494, 48930, 123, 52676, 78870, 53758, 86969, 6745], [10696, 96740, 8521, 97319, 88588, 7634, 36159, 97206, 70536, 16244, 91329]]
#median_solver(S)
fs = median_func(S)


#fs.print_table()
print(fs.optimum(fs.bs,fs.bs))
#2506306.0
#2506397.0  2506397
#print("naive",fs.naive_minimum_fix_a(fs.bs,1))

f = fs.fs[0]

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
