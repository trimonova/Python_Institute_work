class Adder:
    def __init__(self, x):
        self.data = x
    def __add__(self, other):
        return ListAdder(4).add(self.data, other)
    def add(self, y, x):
        print('Not Implemented')

class ListAdder(Adder):
    def add(self, a, b):
        print(a+b)

class DictAdder(Adder):
    def add(self, x, y):
        y.update(x)
        print(y)

e1 = Adder(5)
e2 = ListAdder(6)
e3 = DictAdder(7)

e1.add(5,6)
e2.add([1,2,3], [4,5,6])
e3.add({1:3, 'a': 9}, {'b': 7, 'v': 'a'})

print(e1+6)




