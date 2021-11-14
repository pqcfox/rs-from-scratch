K = 8
MOD = [1, 0, 0, 0, 1, 1, 0, 1, 1]


class QRFiniteField:
    def __init__(self, coeffs):
        self.coeffs = coeffs

    def __str__(self):
        # get list of exponents, then monomials
        exps = range(K-1, -1, -1)
        monos = [f'x^{exp}' if exp != 0 else '1'
                 for coeff, exp in zip(self.coeffs, exps)
                 if coeff == 1]

        # join together with plusses
        poly = ' + '.join(monos)
        return poly

    def __repr__(self):
        # add some context for repr
        return f'<value {self} in field GF(2^{K})>'

    @staticmethod
    def xor(a, b):
        # just xor each element
        return [a_elt ^ b_elt for a_elt, b_elt in zip(a, b)]

    def __add__(self, other):
        # adding is just xoring
        coeffs = xor(self.coeffs, other.coeffs)
        return QRFiniteField(coeffs)

    def __mul__(self, other):
        # create a product with all zeros
        prod = [0 for _ in range(2 * (K-1) + 1)]

        # perform raw multiplication
        exps = range(K-1, -1, -1)

        # go through all self's terms
        for exp, coeff in zip(exps, self.coeffs):
            # if a term is nonzero
            if coeff == 1:

                # create a partial sum by shifting other's coeffs
                part_sum = other.coeffs + exp * [0]
                part_sum = (len(prod) - len(part_sum)) * [0] + part_sum

                # xor it into the product
                prod = self.xor(prod, part_sum)

        # long division to reduce
        while True:
            # find the first nonzero coeff in prod
            prod_exps = range(2 * (K-1), -1, -1)
            prod_degree = max([exp for exp, coeff in zip(prod_exps, prod)
                               if coeff == 1])

            print(prod)
            print(prod_exps)

            # if we're fully reduced, break
            if prod_degree < K:
                break

            # shift our mod to reduce out the top degree
            shift_mod = MOD + [0] * (prod_degree - K)
            shift_mod = (len(prod) - len(shift_mod)) * [0] + shift_mod

            print(shift_mod)
            print('----------')

            # xor it out of the product
            prod = self.xor(prod, shift_mod)

        # chop off the leading bits (now zero) and return
        prod = prod[-K:]
        return QRFiniteField(prod)


if __name__ == '__main__':
    a = QRFiniteField([0, 1, 0, 1, 0, 0, 1, 1])
    b = QRFiniteField([1, 1, 0, 0, 1, 0, 1, 0])

    print(repr(a * b))
