import unittest
import monge
class TestStringMethods(unittest.TestCase):
    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

    def test_stuff(self):
        # test piecewise linear unimodal
        F1 = monge.PiecewiseLinearUnimodal(0.0, [0.0, 1.0, 2.0, 3.0, 4.0], [1, 1, 2, 1, 1], -3)
        F2 = monge.PiecewiseLinearUnimodal(0.0, [0.0, 1.0, 2.0], [1, 3, 1], -2)

        # F1.print()
        assert (F1.slopes() == [-3, -2, -1, 1, 2, 3])
        assert (F1.values() == [0.0, -2.0, -3.0, -2.0, 0.0])

        # F2.print()
        assert (F2.slopes() == [-2, -1, 2, 3])
        assert (F2.values() == [0.0, -1.0, 1.0])
        F1.print()

    def test_median(self):
        S = [[6, 6, 6], [11, 11, 19], [14, 14, 6], [19, 19, 19]]
        assert (monge.median_solver(S)==32.0)

    def test_pure_linear_function(self):
        # a function is of the form ax+b, we solve it by considering it
        a = 3.0
        b = -1
        f = monge.PiecewiseLinearUnimodal(b,[0.0],[0.0],a)
        assert (f.evaluate(3.0) == 8.0)

    def test_flattern(self):
        S = [[6, 6, 6], [11, 11, 19], [14, 14, 6], [19, 19, 19]]
        f = monge.median_func(S)

        g = f.flattern()
        for x in f.bs:
            assert (f.evaluate(x), g.evaluate(x))

if __name__ == '__main__':
    unittest.main()