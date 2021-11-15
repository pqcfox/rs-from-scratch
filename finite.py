K = 8
MOD = [1, 0, 0, 0, 1, 1, 0, 1, 1]


def _degree(poly):
    exps = range(len(poly) - 1, -1, -1)
    degree = max([exp for exp, coeff in zip(exps, poly)
                  if coeff == 1])
    return degree


def _add(poly_a, poly_b):
    # just xor each element
    return [coeff_a ^ coeff_b for coeff_a, coeff_b in zip(poly_a, poly_b)]


def _multiply(poly_a, poly_b):
    # create a product with all zeros
    prod = [0 for _ in range((len(poly_a) - 1) + (len(poly_b) - 1) + 1)]

    # perform raw multiplication
    exps = range(K-1, -1, -1)

    # go through all self's terms
    for exp, coeff_a in zip(exps, poly_a):
        # if a term is nonzero
        if coeff_a == 1:

            # create a partial sum by shifting other's coeffs
            part_sum = poly_b + exp * [0]
            part_sum = (len(prod) - len(part_sum)) * [0] + part_sum

            # add it into the product
            prod = _add(prod, part_sum)

    return prod


def _reduce(poly, mod):
    # get mod degree
    mod_degree = _degree(mod)

    # remove leading zeroes from mod
    mod = mod[-(mod_degree + 1):]

    # create quotient
    quot = [0 for _ in range(len(poly))]

    # long division to reduce
    while True:
        # get poly degree
        poly_degree = _degree(poly)

        # if we're fully reduced, break
        if poly_degree < mod_degree:
            break

        # shift our mod to reduce out the highest monomial
        shift_amount = poly_degree - mod_degree
        shift_mod = MOD + [0] * shift_amount

        # pad out the mod to match poly's length
        shift_mod = (len(poly) - len(shift_mod)) * [0] + shift_mod

        # subtract it out of the product
        poly = _add(poly, shift_mod)

        # add corresponding value to quotient
        quot[-shift_amount] = 1

    return quot, poly


class QRFiniteField:
    def __init__(self, coeffs):
        self.coeffs = coeffs

    def __str__(self):
        # get list of exponents, then monomials
        exps = range(K - 1, -1, -1)
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
        # adding is just adding
        coeffs = _add(self.coeffs, other.coeffs)
        return QRFiniteField(coeffs)

    def __mul__(self, other):
        # do a raw multiply
        prod = _multiply(self.coeffs, other.coeffs)

        # reduce the product
        _, prod = _reduce(prod, MOD)

        # chop off the leading bits (now zero) and return
        prod = prod[-K:]
        return QRFiniteField(prod)


    def inv(self):
        # use extended euclidean method
        # set initial values
        r_last, r_cur = self.coeffs, MOD
        s_last = [0 for _ in range(K)]
        s_last[-1] = 1
        s_cur = [0 for _ in range(K)]

        # iterate until we're done
        while True:
            q_cur, r_next = _reduce(r_last, r_cur)
            if r_next == 0:
                break
            else:
                r_last, r_cur  = r_cur, r_next
                s_last, s_cur = s_cur, _add(s_last, _multiply(q_cur, s_cur))

        # the inverse is stored in s_cur
        return QRFiniteField(s_cur)


if __name__ == '__main__':
    a = QRFiniteField([0, 1, 0, 1, 0, 0, 1, 1])
    b = QRFiniteField([1, 1, 0, 0, 1, 0, 1, 0])

    print(a * b)
