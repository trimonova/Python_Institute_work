class CardHolder:
    acctlen = 8 # Данные класса
    retireage = 59.5
    def __init__(self, acct, name, age, addr):
        self.acct = acct # Данные экземпляра
        self.name = name # Эти операции вызывают методы записи свойств
        self.age = age # К именам вида __X добавляется имя класса
        self.addr = addr # addr – неуправляемый атрибут
        # remain – не имеет фактических данных
    def getName(self):
        return self.__name

    def setName(self, value):
        value = value.lower().replace(' ', '_')
        self.__name = value
    name = property(getName, setName)
    def getAge(self):
        return self.__age
    def setAge(self, value):
        if value < 0 or value > 150:
            raise ValueError('invalid age')
        else:
            self.__age = value
    age = property(getAge, setAge)
    def getAcct(self):
        return self.__acct[:-3] + '***'
    def setAcct(self, value):
        value = value.replace('-', '')
        if len(value) != self.acctlen:
            raise TypeError('invald acct number')
        else:
            self.__acct = value
    acct = property(getAcct, setAcct)
    def remainGet(self): # Можно было бы реализовать как
        return self.retireage - self.age # метод, если нигде не используется
    remain = property(remainGet) # как атрибут

bob = CardHolder('1234-5678', 'Bob Smith', 40, '123 main st')
print(bob.acct, bob.name, bob.age, bob.remain, bob.addr, sep='/')
bob.name = 'Bob Q. Smith'
bob.age = 50
bob.acct = '23-45-67-89'
print(bob.acct, bob.name, bob.age, bob.remain, bob.addr, sep='/')
sue = CardHolder('5678-12-34', 'Sue Jones', 35, '124 main st')
print(sue.acct, sue.name, sue.age, sue.remain, sue.addr, sep= '/')
try:
    sue.age = 200
except:
    print('Bad age for Sue')
try:
    sue.remain = 5
except:
    print('Can’t set sue.remain')
try:
    sue.acct = '1234567'
except:
    print('Bad acct for Sue')