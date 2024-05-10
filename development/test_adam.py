def function(**kwargs):
    a, b, c, d = kwargs['a'], kwargs['b'], kwargs['c'], kwargs['d']
    print(a,b,c,d)

def wrapper(foo):

    def foo_out(**kwargs):

        foo(**kwargs, c=0, d=0)

    return foo_out

out = wrapper(function)

out(a=1, b=2)