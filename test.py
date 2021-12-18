def divideByTwo(val):
    return val / 2


def random():
    a = 10
    if a > 5:
        a = divideByTwo(a)
    return a


a = random()
print(a)
