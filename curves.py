import math

# "Slides" refers to the PDF file "Block 4 Part 2"

# Using a constant variable to change the value for all occurrences centrally.
POINT_AT_INFINITY = "point at infinity"


class CurvesBase:
    """An abstract parent class for the specific curve types containing helper methods and common requirements."""
    coefficient1_name = None
    coefficient2_name = None

    coordinate1_name = "x"
    coordinate2_name = "y"

    def __init__(self, p, coefficient1, coefficient2, x, y) -> None:

        # Checking for possible errors.
        #   Check if all arguments are integers or if the coordinates are the point at infinity otherwise.
        if not isinstance(p, int) or not isinstance(coefficient1, int) or not isinstance(coefficient2, int) or (not isinstance(x, int) and x != POINT_AT_INFINITY) or (not isinstance(y, int) and y != POINT_AT_INFINITY):
            raise TypeError("Domain parameters and coordinates have to be integers.")
        if p <= 3:
            raise ValueError("Prime number p has to be bigger than 3.")

        coefficient1 = coefficient1 % p
        coefficient2 = coefficient2 % p

        if x != POINT_AT_INFINITY and y != POINT_AT_INFINITY:
            x = x % p
            y = y % p

        # Double check: When using modulo p, then this should never be true.
        if coefficient1 not in range(p) or coefficient2 not in range(p):
            raise ValueError("Coefficients have to be within GF(p). Choose an positive integer between 0 and p-1.")

        # Check if p is not prime.
        CurvesBase.if_not_prime_raise_error(p)

        # After checking the values, they should meet all common requirements of the 4 curve types.
        self.p = p
        self.coefficient1 = coefficient1
        self.coefficient2 = coefficient2
        self.x = x
        self.y = y

    @staticmethod
    def if_not_prime_raise_error(p):
        """Raise a ValueError if the given input argument is not a prime number.
        Using Fermat's Little Theorem: If a**n mod n = a, then n is likely prime. If this is not true, n is not prime.
        E.g.: p=2**251-1 isn't prime, but passes pow(2, p, p)==2, yet we can test with more bases as needed.
        """
        if pow(2, p, p) != 2 or pow(3, p, p) != 3:
            raise ValueError("Given value for p is not a prime number. Please choose a prime number bigger than 3. Given p = " + str(p))
        # else: pow returns 2, so p could be prime (but not guaranteed).

    @staticmethod
    def get_legendre_symbol(a, p):
        """Calculate and return legendre symbol of given argument in given prime (field)."""
        a = a % p
        if a % p == 0:
            return 0
        legendre_symbol = pow(a, (p - 1) // 2, p)
        # Results are {0, 1, −1} in Fp, -1 is congruent to p-1.
        if legendre_symbol == p - 1:  # Don't omit for returning p-1, because other functions test for -1.
            legendre_symbol = -1
        return legendre_symbol

    @staticmethod
    def is_quadratic_residue(a, p):
        """Call get_legendre_symbol(), interpret and return if a is a qr in p."""
        legendre_symbol = CurvesBase.get_legendre_symbol(a, p)
        # Following https://en.wikipedia.org/wiki/Legendre_symbol#Definition last accessed on 2022-07-26
        if legendre_symbol == 1:
            return True
        elif legendre_symbol == -1:
            return False
        elif legendre_symbol == 0:
            return None
        else:
            raise ValueError("Check if p is prime, otherwise problem in is_quadratic_residue() or legendre_symbol().")

    @staticmethod
    def get_modular_multiplicative_inverse(factor, p):
        """Calculate multiplicative inverse of given factor and given prime.
        Uses Fermat's little theorem.
        p has to be prime. 
        Should be the case when invoked from object, because of the check in init of CurveBase.
        Control the calculated value for the inverse with assert.
        """
        if math.gcd(factor, p) != 1:
            raise ValueError(f"The factor {factor} and p {p} have to be coprime integers.")
        # Throw error if p not prime.
        CurvesBase.if_not_prime_raise_error(p)
        # using Fermat's little theorem:
        inverse = pow(factor, p - 2, p)
        # Check:
        assert factor * inverse % p == 1
        return inverse

    def print_transformation_results(self):
        """Generic method for all specific curves."""
        print("\t Results:")
        output_dic = {
            self.coefficient1_name: self.coefficient1,
            self.coefficient2_name: self.coefficient2,
            self.coordinate1_name: self.x,
            self.coordinate2_name: self.y
        }
        # Print out values.
        for key in output_dic:
            print(f"The value for {key} = {output_dic[key]}")

        # Print out the equation of resulting curve.
        print(self.equation)

    def get_modular_sqrt(self, quadrat):
        """Calculate and return modular square root using Tonelli–Shanks algorithm.
        https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm for information (including core ideas and pseudocode) last accessed on 2022-07-13.
        Faster than checking/brute-forcing values by squaring v mod p and check if it matches quadrat.
        """

        # Check if quadrat is not a quadratic residue.
        if self.get_legendre_symbol(quadrat, self.p) != 1:
            raise ValueError(f"Parameter {quadrat} is not a qr, so one cannot find a sqrt.")
        elif quadrat == 0:
            return 0
        elif self.p % 4 == 3:  # simplified Tonelli–Shanks https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm#The_algorithm last accessed on 2022-07-13, found in the last line
            return pow(quadrat, ((self.p + 1) // 4), self.p)
        elif self.p % 8 == 5:  # https://www.rieselprime.de/ziki/Modular_square_root last accessed on 2022-07-14
            v = pow(2 * quadrat, (self.p - 5) // 8, self.p)
            i = 2 * quadrat * v ** 2 % self.p
            root = quadrat * v * (i - 1) % self.p
            return root

        even_number = self.p - 1  # Is always a composite.
        exponent = 0  # Exponent for 2**exponent.
        while even_number % 2 == 0:
            even_number = even_number // 2
            exponent = exponent + 1
        # At the end of the while loop, value of even_number is odd.
        odd_number = even_number

        # Find a qnr z.
        z = 2
        while self.get_legendre_symbol(z, self.p) != -1:
            z += 1
        # Now self.get_legendre_symbol(z , self.p) is -1.
        quadratic_non_residue = z

        M = exponent
        R = quadrat ** ((odd_number + 1) // 2) % self.p
        c = quadratic_non_residue ** odd_number % self.p
        t = quadrat ** odd_number % self.p
        while True:
            if t == 0:
                return 0
            elif t == 1:
                return int(R)

            i = 0

            for i in range(M):
                t = t ** (2 ** i) % self.p
                if t == 1:
                    break

            b = pow(c, (2 ** (M - i - 1)), self.p)
            M = i
            c = b ** 2 % self.p
            t = t * b ** 2 % self.p
            R = R * b % self.p


class ShortWeierstrass(CurvesBase):
    """A class for Short Weierstrass curves with init, transformations and helper methods."""
    curve_name = "Short Weierstrass"
    coefficient1_name = "a"
    coefficient2_name = "b"

    def __init__(self, p, coefficient1, coefficient2, x, y) -> None:
        super().__init__(p, coefficient1, coefficient2, x, y)
        self.a = self.coefficient1
        self.b = self.coefficient2

        self.available_transformations = {
            "to m": self.trans_sw_to_m,
            "to ted": self.trans_sw_to_ted,
            "to ed": self.trans_sw_to_ed
        }

        # Checking the coefficients.
        if (4 * self.a ** 3 + 27 * self.b ** 2) % self.p == 0:
            # In this case the discriminant dissolves.
            raise ValueError("Coefficients provided result in an singular curve. The curve requires to be non-singular. Please review the coefficients.")
            # Maybe paste the equation of the if statement in the error message. 

        # equation: y**2 = x**3 + a * x + b
        self.equation = f"Equation: y**2 = x**3 + {self.a} * x + {self.b}"

    def get_zero_of_sw_function(self):
        # Find a zero of a sw function zof.
        for x in range(self.p):
            RHS = (x ** 3 + self.a * x + self.b) % self.p
            if RHS == 0:
                return x
        raise ValueError("Could not find a zero of a function for the given domain parameters.")

    def trans_sw_to_m(self, print_result=False):
        # Slides page 9.
        # 2.1 in meeting notes.
        # Check conditions for the transformation:
        #   1. x**3 + a*x + b must have one zero of a function zof.
        zof = self.get_zero_of_sw_function()

        #   2. 3 * zof**2 + a is a qr
        check_number = (3 * zof ** 2 + self.a) % self.p
        if not self.is_quadratic_residue(check_number, self.p):
            raise ValueError("Condition for sw to m transformation failed. 3 * (the smallest existing zero of a function) **2 + a has to be a qr in Fp.")

        # Get s for further transformation:
        s_inverse = self.get_modular_sqrt(check_number)
        s = CurvesBase.get_modular_multiplicative_inverse(s_inverse, self.p)

        # Coefficients:
        A = 3 * zof * s % self.p
        B = s % self.p

        # Points:
        u = s * (self.x - zof) % self.p
        v = s * self.y % self.p

        result_curve = Montgomery(self.p, A, B, u, v)
        if print_result:
            result_curve.print_transformation_results()

        return self.p, A, B, u, v

    def trans_sw_to_ted(self, print_result=False):
        # 2.3 in meeting notes.
        # Transform sw to m and m to ted.
        # Not explicitly given in slides.

        mont_params = self.trans_sw_to_m()
        intermediate_curve = Montgomery(*mont_params)

        result_params = intermediate_curve.trans_m_to_ted(print_result)

        return result_params

    def trans_sw_to_ed(self, print_result=False):
        # 2.2 in meeting notes.
        # Transform sw to m and m to ted and ted to ed.
        # Not explicitly given in slides.

        mont_params = self.trans_sw_to_m()
        intermediate_curve_m = Montgomery(*mont_params)

        ted_params = intermediate_curve_m.trans_m_to_ted()
        intermediate_curve_ted = TwistedEdwards(*ted_params)

        result_params = intermediate_curve_ted.trans_ted_to_ed(print_result)

        return result_params


class Montgomery(CurvesBase):
    """A class for Montgomery curves with init and transformations."""
    curve_name = "Montgomery"
    coefficient1_name = "A"
    coefficient2_name = "B"

    coordinate1_name = "u"
    coordinate2_name = "v"

    def __init__(self, p, coefficient1, coefficient2, x, y) -> None:
        super().__init__(p, coefficient1, coefficient2, x, y)
        self.A = self.coefficient1
        self.B = self.coefficient2

        self.available_transformations = {
            "to ted": self.trans_m_to_ted,
            "to sw": self.trans_m_to_sw,
            "to ed": self.trans_m_to_ed
        }

        # Checking the coefficients.
        if self.A == 2 or self.A == -2:
            raise ValueError("Coefficient A cannot be 2 or -2.")
        if self.B == 0:
            raise ValueError("Coefficient B cannot be 0.")
        if self.B * (self.A ** 2 - 4) % self.p == 0:
            raise ValueError("At least one of the coefficients has to change.")

        # equation: B * y**2 = (x**3 + A*x**2 + x) % p
        self.equation = f"Equation: {self.B} * v**2 = u**3 + {self.A} * u**2 + u"

    def trans_m_to_ted(self, print_result=False):
        # Slides page 4.
        # 3.2 in meeting notes.

        # Check for isomorphism.
        #   For a, which has to be a qr:
        no_isomorphism_error_message = "Isomorphism, the condition for this transformation, is not given. Please review A and B."
        if not CurvesBase.is_quadratic_residue((self.A + 2) * self.get_modular_multiplicative_inverse(self.B, self.p), self.p):
            raise ValueError(no_isomorphism_error_message + " Result of (A + 2)/B isn't a quadratic residue.")

        #   For d, which has to be a non-qr:
        if CurvesBase.is_quadratic_residue((self.A - 2) * self.get_modular_multiplicative_inverse(self.B, self.p), self.p):
            raise ValueError(no_isomorphism_error_message + " Error: Result of (A - 2)/B is a quadratic residue.")

        # Montgomery naming convention for coordinates: u and v for TEd x and y
        u = self.x
        v = self.y

        # Transformation
        # Points
        if u == POINT_AT_INFINITY and v == POINT_AT_INFINITY:
            x = 0
            y = 1
        elif u == 0 and v == 0:
            x = 0
            y = -1 % self.p
        else:
            x = u * CurvesBase.get_modular_multiplicative_inverse(v, self.p) % self.p
            y = (u - 1) * CurvesBase.get_modular_multiplicative_inverse(u + 1, self.p) % self.p

        # Coefficients:
        a = (self.A + 2) * CurvesBase.get_modular_multiplicative_inverse(self.B, self.p) % self.p
        d = (self.A - 2) * CurvesBase.get_modular_multiplicative_inverse(self.B, self.p) % self.p

        result_curve = TwistedEdwards(self.p, a, d, x, y)
        if print_result:
            result_curve.print_transformation_results()

        return self.p, a, d, x, y

    def trans_m_to_sw(self, print_result=False):
        # Slides page 7.
        # 1.1 in meeting notes.
        # Works always, so no need to check for conditions.

        # Points:
        # Montgomery naming convention for coordinates: u and v to SW x and y
        u = self.x
        v = self.y

        x = (u * CurvesBase.get_modular_multiplicative_inverse(self.B, self.p) + self.A * CurvesBase.get_modular_multiplicative_inverse(3 * self.B, self.p)) % self.p
        y = v * CurvesBase.get_modular_multiplicative_inverse(self.B, self.p) % self.p

        # Coefficients: 
        a = (3 - self.A ** 2) * CurvesBase.get_modular_multiplicative_inverse(3 * self.B ** 2, self.p) % self.p
        b = (2 * self.A ** 3 - 9 * self.A) * CurvesBase.get_modular_multiplicative_inverse(27 * self.B ** 3, self.p) % self.p

        result_curve = ShortWeierstrass(self.p, a, b, x, y)
        if print_result:
            result_curve.print_transformation_results()

        return self.p, a, b, x, y

    def trans_m_to_ed(self, print_result=False):
        # 3.1 in meeting notes.
        # Transform m to ted and ted to ed.
        # Not explicitly given in slides.

        ted_params = self.trans_m_to_ted()
        intermediate_curve_ted = TwistedEdwards(*ted_params)

        result_params = intermediate_curve_ted.trans_ted_to_ed(print_result)

        return result_params


class Edwards(CurvesBase):
    """A class for Edwards curves with init and transformations."""
    curve_name = "Edwards"
    coefficient1_name = "c"
    coefficient2_name = "d"

    def __init__(self, p, coefficient1, coefficient2, x, y) -> None:
        super().__init__(p, coefficient1, coefficient2, x, y)
        self.c = self.coefficient1
        self.d = self.coefficient2

        self.available_transformations = {
            "to m": self.trans_ed_to_m,
            "to sw": self.trans_ed_to_sw,
            "to ted": self.trans_ed_to_ted
        }

        # Check the coefficient d.
        if self.d == 0 or self.d == 1:
            raise ValueError("Coefficient d cannot be 0 or 1.")
        # Coefficient d has to be a qnr for the point addition to be complete. https://safecurves.cr.yp.to/equation.html last accessed on 2022-07-21
        if CurvesBase.is_quadratic_residue(self.d, self.p):
            raise ValueError("Coefficient d has to be a qnr in Fp. The chosen d is a qr in Fp. Please adjust d or p.")
        # d(1-d) has to be nonzero in F_p. https://safecurves.cr.yp.to/equation.html last accessed on 2022-07-21
        if self.d * (1 - self.d) % self.p == 0:
            raise ValueError("Please review the given d or p. Note that d(1-d) has to be nonzero mod p.")

        # Check the coefficients c and d.
        # cd(1 − dc**4) cannot be 0. https://eprint.iacr.org/2007/286.pdf page 2, last accessed on 2022-07-17
        if (self.c * self.d * (1 - self.c ** 4 * self.d)) % self.p == 0:  # also on https://en.wikipedia.org/wiki/Edwards_curve#Definition last accessed on 2022-07-15
            raise ValueError("Coefficients c and d can't have these exact values at the same time.")

        # equation: x**2 + y**2 = c**2 * (1 + d * x**2 * y**2)
        self.equation = f"Equation: x**2 + y**2 = {self.c}**2 * (1 + {self.d} * x**2 * y**2)"

    def trans_ed_to_m(self, print_result=False):
        # Not explicitly given in slides.
        # 4.1 in meeting notes.
        # https://safecurves.cr.yp.to/verify.html last accessed on 2022-07-21
        # Requires that c = 1.
        if self.c != 1:
            raise ValueError("Coefficient c should be 1 for this transformation.")

        A = 2 * (1 + self.d) * CurvesBase.get_modular_multiplicative_inverse(1 - self.d, self.p) % self.p
        B = 4 * CurvesBase.get_modular_multiplicative_inverse(1 - self.d, self.p) % self.p

        u = (1 + self.y) * CurvesBase.get_modular_multiplicative_inverse(1 - self.y, self.p) % self.p
        v = ((1 + self.y) * CurvesBase.get_modular_multiplicative_inverse(1 - self.y, self.p)) * CurvesBase.get_modular_multiplicative_inverse(self.x, self.p) % self.p

        resulting_curve = Montgomery(self.p, A, B, u, v)
        if print_result:
            resulting_curve.print_transformation_results()

        return self.p, A, B, u, v

    def trans_ed_to_sw(self, print_result=False):
        # 1.2 in meeting notes.
        # Transform ed to m and m to sw.
        # Not explicitly given in slides.

        mont_params = self.trans_ed_to_m()
        intermediate_curve = Montgomery(*mont_params)

        result_params = intermediate_curve.trans_m_to_sw(print_result)

        return result_params

    def trans_ed_to_ted(self, print_result=False):
        # Not explicitly given in slides.
        # Given in "Block 4 Part 1" page 8.
        # Not in meeting notes.

        # Requires that c = 1.
        if self.c != 1:
            raise ValueError("Coefficient c should be 1 for this transformation.")

        # Find a value for coefficient a, so that coefficient a is a qr in Fp.
        a = None
        for i in range(2, self.p):
            if CurvesBase.is_quadratic_residue(i, self.p):
                a = i
                break
        if a is None:
            raise ValueError("Could not find a qr in Fp.")

        # self.var are parameters from the Edwards curve.
        ted_d = self.d * a % self.p

        # Points:
        #   ted_x = x/sqrt(a)
        ted_x = self.x * self.get_modular_multiplicative_inverse(self.get_modular_sqrt(a), self.p) % self.p
        ted_y = self.y

        result_curve = TwistedEdwards(self.p, a, ted_d, ted_x, ted_y)
        if print_result:
            result_curve.print_transformation_results()

        return self.p, a, ted_d, ted_x, ted_y


class TwistedEdwards(CurvesBase):
    """A class for Twisted Edwards curves with init and transformations."""
    curve_name = "Twisted Edwards"
    coefficient1_name = "a"
    coefficient2_name = "d"

    def __init__(self, p, coefficient1, coefficient2, x, y) -> None:
        super().__init__(p, coefficient1, coefficient2, x, y)
        self.a = self.coefficient1
        self.d = self.coefficient2

        # https://eprint.iacr.org/2008/013.pdf § Definition 2.1 last accessed on 2022-07-27
        # Check that a and d are nonzero.
        if self.a == 0 or self.d == 0:  # also https://en.wikipedia.org/wiki/Twisted_Edwards_curve#Definition last accessed on 2022-07-20
            raise ValueError("Coefficients a and d should be nonzero elements (of Fp).")

        # Check that a and d are distinct.
        # https://eprint.iacr.org/2008/013.pdf § Definition 2.1 last accessed on 2022-07-27
        # also https://en.wikipedia.org/wiki/Twisted_Edwards_curve#Definition last accessed on 2022-07-26
        if self.a == self.d:
            raise ValueError("Coefficients a and d should be distinct elements (of Fp).")

        self.available_transformations = {
            "to m": self.trans_ted_to_m,
            "to sw": self.trans_ted_to_sw,
            "to ed": self.trans_ted_to_ed
        }

        # Possibly the same limit for d not being 0 or 1. 

        # equation: a * x**2 + y**2 = 1 + d * x**2 * y**2
        self.equation = f"Equation: {self.a} * x**2 + y**2 = 1 + {self.d} * x**2 * y**2"

    def trans_ted_to_m(self, print_result=False):
        # 5.1 in meeting notes.
        # Page 3 in Slides.
        # Also https://de.wikipedia.org/wiki/Edwards-Kurve last accessed on 2022-07-22.
        # Input arguments can be replaced in the func with self.var to make function work with less input. 
        # Transform points and domain parameters to a Montgomery curve.

        # Addition has to be complete.         
        # Is complete if a is a qr and d is not a qr in Fp.
        if not CurvesBase.is_quadratic_residue(self.a, self.p):
            raise ValueError("Coefficient a has to be a quadratic residue in Fp.")
        if CurvesBase.is_quadratic_residue(self.d, self.p):
            raise ValueError("Coefficient d can't be a quadratic residue in Fp. Currently d = " + str(self.d))

        # Points:
        # u is x Coordinate on montgomery, v ist y
        if self.x == 0 and self.y == 1:
            u = POINT_AT_INFINITY
            v = POINT_AT_INFINITY
        elif self.x == 0 and (self.y == -1 or self.y == self.p - 1):  # Could be -1 mod p so also p-1.
            u = 0
            v = 0
        else:
            # u = (1+y)/(1-y)
            u = (1 + self.y) * CurvesBase.get_modular_multiplicative_inverse(1 - self.y, self.p) % self.p

            # v = (1+y)/(x*(1-y))
            v = (1 + self.y) * CurvesBase.get_modular_multiplicative_inverse(self.x * (1 - self.y), self.p) % self.p

        # Coefficients:
        # A = 2 * (a+d)/(a-d)
        A = 2 * (self.a + self.d) * CurvesBase.get_modular_multiplicative_inverse(self.a - self.d, self.p)
        A = A % self.p

        # B = 4/(a-d)
        B = 4 * CurvesBase.get_modular_multiplicative_inverse(self.a - self.d, self.p)
        B = B % self.p

        resulting_curve = Montgomery(self.p, A, B, u, v)
        if print_result:
            resulting_curve.print_transformation_results()

        return self.p, A, B, u, v

    def trans_ted_to_sw(self, print_result=False):
        # 1.3 in meeting notes.
        # Transform ted to m and m to sw.
        # Not explicitly given in slides.

        mont_params = self.trans_ted_to_m()
        intermediate_curve = Montgomery(*mont_params)

        result_params = intermediate_curve.trans_m_to_sw(print_result)

        return result_params

    def trans_ted_to_ed(self, print_result=False):
        # Not explicitly given in slides.
        # Inverted method of given one in "Block 4 Part 1" page 8.
        # Not in meeting notes.

        # Sets c = 1.
        c = 1
        ed_d = self.d * CurvesBase.get_modular_multiplicative_inverse(self.a, self.p) % self.p

        # get_modular_sqrt checks if self.a is not a quadratic residue.
        ed_x = self.x * self.get_modular_sqrt(self.a) % self.p
        ed_y = self.y

        resulting_curve = Edwards(self.p, c, ed_d, ed_x, ed_y)
        if print_result:
            resulting_curve.print_transformation_results()

        return self.p, c, ed_d, ed_x, ed_y
