import time
def timer(label='', trace=True): # Аргументы декоратора: сохраняются
    class Timer:
        def __init__(self, func): # На этапе декорирования сохраняется
            self.func = func # декорируемая функция
            self.alltime = 0
        def __call__(self, *args, **kargs): # При вызове: вызывается оригинал
            start = time.clock()
            result = self.func(*args, **kargs)
            elapsed = time.clock() - start
            self.alltime += elapsed
            if trace:
                format = '%s %s: %.5f, %.5f'
                values = (label, self.func.__name__, elapsed, self.alltime)
                print(format % values)
            return result
    return Timer

@timer(label='[CCC]==>')
def listcomp(N): # То же, что и listcomp = timer(...)(listcomp)
    return [x * 2 for x in range(N)] # listcomp(...) вызовет Timer.__call__

@timer(trace=True, label='[MMM]==>')
def mapcall(N):
    return list(map((lambda x: x * 2), range(N)))
for func in (listcomp, mapcall):
    print('')

    result = func(5) # Хронометраж вызова, всех вызовов, возвращаемое значение
    func(50000)
    func(500000)
    func(1000000)
    print(result)
    print('allTime = %s' % func.alltime) # Общее время всех вызовов
print('map/comp = %s' % round(mapcall.alltime / listcomp.alltime, 3))