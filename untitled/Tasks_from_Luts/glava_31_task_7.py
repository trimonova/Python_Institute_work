class Lunch:
    def __init__(self):
        self.client = Customer()
        self.oficiant = Employee()
        self.zakaz = Food()
    def order(self, foodName):
        self.client.placeOrder(foodName, self.oficiant, self.zakaz)
    def result(self):
        self.client.printFood(self.zakaz.menu)


class Customer:
    #def __init__(self):
    def placeOrder(self, foodName, employee, zakaz):
        employee.takeOrder(foodName, zakaz)
    def printFood(self, zakaz):
        print(zakaz)


class Employee:
    def takeOrder(self, foodName, zakaz):
        zakaz.whole_menu(foodName)


class Food:

    def __init__(self):
        self.menu = []
    def whole_menu(self, name):
        self.menu.append(name)


lunch = Lunch()
lunch.order('pizza')
lunch.order('burito')
lunch.result()


lunch2 = Lunch()
lunch2.order('pasta')
lunch2.order('gril')
lunch2.result()
