class ListInstance:
    def __str__(self):
        return '<Instance of {}, address {}:\n{}>'.format(self.__class__.__name__, id(self), self.__attrnames())
    def __attrnames(self):
        result = ''
        for attr in sorted(self.__class__.__bases__): # Словарь атрибутов
            result += '\tname {}={}\n'.format(attr, attr.__name__ )
        return result

class Spam(ListInstance): # Наследует метод __str__
    def __init__(self):
        self.data1 = 'food'

x = Spam()
print(x)