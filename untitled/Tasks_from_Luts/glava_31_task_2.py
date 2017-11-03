class MyList:
    def __init__(self, list_1 = []):
        self.list_1 = list_1
    def __add__(self, other):
        return self.list_1 + other
    def __getitem__(self, item):
        return self.list_1[item]
    #def __iter__(self):
        #return next(MyList)
    def append (self, element):
        self.list_1.append(element)
        return self.list_1
    def sort (self):
        return self.list_1.sort()
    #def __str__(self):
        #return 'MyList {}'.format(self.list_1)

#a = MyList()
#b = a + [5,8,7]
#print(b[0:2])
#
#print(b[1])
#
#for i in b:
#    print(i)
#
#b.append(11)
#print(b)
#b.sort()
#print(b)

class MylistSub(MyList):
    count_add = 0
    count_getitem = 0
    def __init__(self, list_2 = []):

        self.list_2 = list_2
        MyList.__init__(self, list_2)
    def __add__(self, other):
        MylistSub.count_add += 1
        print('Add is done')
        return MylistSub(MyList.__add__(self, other))

    def __getitem__(self, item):
        MylistSub.count_getitem += 1
        print('Get_item is done')
        return MylistSub(MyList.__getitem__(self, item))
    def __str__(self):
        return 'count_add = {}, count_getitem = {}'.format(self.count_add, self.count_getitem)

c = MyList([1,5,2])
print(c+[2,3])

a = MylistSub([1,5,2])
b = a + [3, 7]
b[1:4]
b = b + [6,8]
a = b[1:3]

print(a)
print(b)



