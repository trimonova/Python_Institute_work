class Meta(BaseException):
    count_person = 0
    def __init__(self, name = '', phone = '', age = ''):
        self.name = name
        self.phone = phone
        self.age = age
        Meta.count_person += 1
    def object(self, person_dict = {}):
        person_dict[self.name] = [self.phone, self.age]
        return person_dict

    def __str__(self):
        return '<{}>'.format(self.object())

    def __getattr__(self, item):
        if item in self.__dict__:
            print('all good')
        else:
            print('all bad')
            raise 'AttributeError <{} is not Meta attribute>'.format(item)

    def __setattr__(self, key, value):
        print(key, value)

a = Meta('Masha', '41-94-22', '25')
a.k = 5



