def divideByTwo(val):
    return val / 2


def random():
    a = 10
    if a > 5:
        a = divideByTwo(a)
    return a


a = random()
print(a)

switch = 0
if seed_direction == 'Forward':
    seed_direction = 'Reverse'
    direction_x = 'Upstream'
    CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
else:
    break_ = True
