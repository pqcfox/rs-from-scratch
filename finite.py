# K is the exponent in our field GL(2^K)
K = 8

# MOD is the coefficients of the polynomial modulus we choose for GL(2^K),
# which here is x^8 + x^4 + x^3 + x + 1
MOD = [1, 0, 0, 0, 1, 1, 0, 1, 1]

# trim leading zeros from a polynomial (list of coefficents, starting with
# the largest-degree monomial)
def _trim(poly):
    if _is_zero(poly):
        return []

    return poly[-(_degree(poly) + 1):]


# test if a polynomial is zero
def _is_zero(poly):
    return all([coeff == 0 for coeff in poly])


# get the degree of a polynomial 
def _degree(poly):
    if _is_zero(poly):
        return None

    # get the highest-degree monomial with non-zero coefficient
    exps = range(len(poly) - 1, -1, -1)
    degree = max([exp for exp, coeff in zip(exps, poly)
                  if coeff == 1])
    return degree


# add two binary polynomials (i.e. as elts of Z/2Z[X])
def _add(poly_a, poly_b):
    # zero-extend the shorter
    len_diff = len(poly_a) - len(poly_b)
    if len_diff > 0:
        poly_b = [0] * len_diff + poly_b
    elif len_diff < 0:
        poly_a = [0] * (-len_diff) + poly_a

    # just xor each coefficient
    return [coeff_a ^ coeff_b for coeff_a, coeff_b in zip(poly_a, poly_b)]


# multiply two binary polynomials
def _multiply(poly_a, poly_b):
    # create a product with all zeros
    prod = [0 for _ in range((len(poly_a) - 1) + (len(poly_b) - 1) + 1)]

    # go through all self's terms
    exps = range(len(poly_a)-1, -1, -1)
    for exp, coeff_a in zip(exps, poly_a):
        # if a term is nonzero
        if coeff_a == 1:

            # create a partial sum by shifting other's coefficients
            part_sum = poly_b + exp * [0]
            part_sum = (len(prod) - len(part_sum)) * [0] + part_sum

            # add it into the product
            prod = _add(prod, part_sum)

    return prod


def _reduce(poly, mod):
    # get mod degree
    mod_degree = _degree(mod)

    # remove leading zeroes from mod
    mod = _trim(mod)

    # create quotient
    quot = [0 for _ in range(len(poly))]

    # long division to reduce
    while True:
        # get poly degree
        poly_degree = _degree(poly)

        # if we're fully reduced, break
        if poly_degree is None or poly_degree < mod_degree:
            break

        # shift mod to cancel out the highest monomial
        shift_amount = poly_degree - mod_degree
        shift_mod = mod + [0] * shift_amount

        # left pad the mod to match poly's length
        shift_mod = (len(poly) - len(shift_mod)) * [0] + shift_mod

        # subtract result out of poly
        poly = _add(poly, shift_mod)

        # add corresponding value to quotient
        quot[-(shift_amount + 1)] = 1

    return quot, poly


class QRFiniteField:
    """
    An element of the finite field GL(2^8).

    Attributes
    ----------
    coeffs : list[int]
        A list of binary coefficients, starting with highest-degree
        monomial, of the element of GL(2^8).
    """

    def __init__(self, coeffs):
        self.coeffs = _trim(coeffs)

    def __str__(self):
        if _is_zero(self.coeffs):
            poly = '0'
        else:
            # get list of exponents, then monomials
            exps = range(len(self.coeffs) - 1, -1, -1)
            monos = [f'x^{exp}' if exp != 0 else '1'
                     for coeff, exp in zip(self.coeffs, exps)
                     if coeff == 1]

            # join together with plusses
            poly = ' + '.join(monos)
        return poly

    def __repr__(self):
        # add some context for repr
        return f'<value {self} in field GF(2^{K})>'

    def __add__(self, other):
        # adding is just adding binary polynomials
        coeffs = _add(self.coeffs, other.coeffs)
        return QRFiniteField(coeffs)

    def __mul__(self, other):
        # multiply as regular binary polynomials
        prod = _multiply(self.coeffs, other.coeffs)

        # reduce the product
        _, prod = _reduce(prod, MOD)

        # return the result
        return QRFiniteField(prod)


    def inv(self):
        """
        Compute the inverse of this element in GL(2^8).

        Returns
        -------
        QRFiniteField:
            The inverse of this element.
        """
        # use extended euclidean method
        # set initial values
        r_last, r_cur = self.coeffs, MOD
        s_last, s_cur = [1], []  # polynomials 1, 0

        # iterate until we're done
        while True:
            # calculate next values
            q_cur, r_next = _reduce(r_last, r_cur)
            s_next = _add(s_last, _multiply(q_cur, s_cur))

            if _is_zero(r_next):
                # quit if no remainder
                break
            else:
                # otherwise, step forward
                r_last, r_cur  = r_cur, r_next
                s_last, s_cur = s_cur, s_next

        # the inverse is stored in s_cur at the end
        return QRFiniteField(s_cur)


if __name__ == '__main__':
    # run a basic test from Wikipedia, Finite Fields
    a = QRFiniteField([0, 1, 0, 1, 0, 0, 1, 1])
    b = QRFiniteField([1, 1, 0, 0, 1, 0, 1, 0])

    assert a.inv().coeffs == b.coeffs
    assert (a * a.inv()).coeffs == [1]
