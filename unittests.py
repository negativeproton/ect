import unittest

from curves import *


# syntax for assertRaises:
# self.assertRaises(ValueError, CurvesBase.if_not_prime_raise_error, 190)
# The argument has to be a value that should trigger the specified error.
# Note how there are no parenthesis and arguments directly after the method, instead being given in the same order after a comma.

# One can test for false negatives with correct values, e.g.
# self.assertRaises(ValueError,CurvesBase.if_not_prime_raise_error, 191)

class TestCurvesMethods(unittest.TestCase):

    def test_if_not_prime_raise_error(self):
        self.assertRaises(ValueError, CurvesBase.if_not_prime_raise_error, -1)
        self.assertRaises(ValueError, CurvesBase.if_not_prime_raise_error, -599999)
        self.assertRaises(ValueError, CurvesBase.if_not_prime_raise_error, 1000000)
        self.assertRaises(TypeError, CurvesBase.if_not_prime_raise_error, 0.0000001)
        self.assertRaises(ValueError, CurvesBase.if_not_prime_raise_error, 190)
        self.assertRaises(TypeError, CurvesBase.if_not_prime_raise_error, "string")
        self.assertRaises(TypeError, CurvesBase.if_not_prime_raise_error, b"05")

        # What should not throw an error:
        CurvesBase.if_not_prime_raise_error(191)

    def test_get_legendre_symbol(self):
        # Case: Result is 0:
        self.assertEqual(CurvesBase.get_legendre_symbol(191, 191), 0)
        self.assertEqual(CurvesBase.get_legendre_symbol(0, 191), 0)

    def test_get_modular_sqrt(self):
        test_curve = CurvesBase(37, 45, 54, 45, 54)
        self.assertEqual(test_curve.get_modular_sqrt(30), 20)

        test_curve = TwistedEdwards(41, 15, 16, 0, -1)
        self.assertEqual(test_curve.get_modular_sqrt(5), 28)

    def test_trans_m_to_ted(self):
        # for Curve25519 (M)
        # Reverting values from https://ed25519.cr.yp.to/ed25519-20110926.pdf § "Choice of curve" and https://math.stackexchange.com/questions/1392277/point-conversion-between-twisted-edwards-and-montgomery-curves both last accessed on 2022-07-24.
        p = 2 ** 255 - 19
        A = 486662
        B = 57896044618658097711785492504343953926634992332820282019728792003956564333285
        u = 9
        v = 46155036877857898950720737868668298259344786430663990124372813544693780678454
        Curve25519 = Montgomery(p, A, B, u, v)

        # to Ed25519 (TEd):
        a = -1 % p
        d = -121665 * 37095705934669439343138083508754565189542113879843219016388785533085940283556 % p
        x = 15112221349535400772501151409588531511454012693041857206046113283949847762202
        y = 46316835694926478169428394003475163141307993866256225615783033603165251855960
        self.assertEqual(Curve25519.trans_m_to_ted(), (p, a, d, x, y))

    def test_trans_m_to_ed(self):
        # https://ed25519.cr.yp.to/ed25519-20110926.pdf § "Choice of curve" last accessed on 2022-07-24.
        # from Curve25519 (M)
        p = 2 ** 255 - 19
        A = 486662
        B = 1
        u = 9
        v = 9
        Curve25519 = Montgomery(p, A, B, u, v)

        # to Edwards:
        c = 1
        d = 20800338683988658368647408995589388737092878452977063003340006470870624536394  # = 21665 * inverse(121666) % p
        x = 9094040566125962849133224048217411091405536248825867518642941381412595940312
        y = 46316835694926478169428394003475163141307993866256225615783033603165251855960
        self.assertEqual(Curve25519.trans_m_to_ed(), (p, c, d, x, y))

    def test_trans_ted_to_m(self):
        test_curve = TwistedEdwards(191, 15, 61, 0, 1)
        self.assertEqual(test_curve.trans_ted_to_m(), (191, 5, 166, POINT_AT_INFINITY, POINT_AT_INFINITY))

        # for Ed25519 (TEd):
        # https://ed25519.cr.yp.to/ed25519-20110926.pdf § "Choice of curve" and https://math.stackexchange.com/questions/1392277/point-conversion-between-twisted-edwards-and-montgomery-curves both last accessed on 2022-07-24.
        p = 2 ** 255 - 19
        a = -1
        d = -121665 * 37095705934669439343138083508754565189542113879843219016388785533085940283556
        x = 15112221349535400772501151409588531511454012693041857206046113283949847762202
        y = 46316835694926478169428394003475163141307993866256225615783033603165251855960
        Ed25519 = TwistedEdwards(p, a, d, x, y)

        # to Curve25519 (M)
        A = 486662
        B = 57896044618658097711785492504343953926634992332820282019728792003956564333285
        u = 9
        v = 46155036877857898950720737868668298259344786430663990124372813544693780678454
        self.assertEqual(Ed25519.trans_ted_to_m(), (p, A, B, u, v))

        # Check init method.
        with self.assertRaises(ValueError):
            my_curve = TwistedEdwards(41, 15, 15, 0, -1)


if __name__ == '__main__':
    unittest.main()

