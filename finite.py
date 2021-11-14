class FiniteField:
    def __init__(self, p, k):
        self.p = p
        self.k = k
        self.q = p**k

    def __call__(self, val):
        return FiniteFieldElt(self, val)

    def __str__(self):
        return f'GF({self.q})'

    def __repr__(self):
        return f'FiniteField({self.p}, {self.k})'


class FiniteFieldElt:
    def __init__(self, val, field):
        self.val = val % field.q
        self.field = field

    def __str__(self):
        return f'{self.val}'

    def __repr__(self):
        return f'<value {self.val} in {self.field}>'


if __name__ == '__main__':
    k = FiniteField(2, 3)
    a = FiniteFieldElt(3, k)
    print(repr(a))
