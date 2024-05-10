'''argument = dict(
    b=2,
    c=0,
    d=4
)

def function(**kwargs):
    lis = []
    for i in kwargs:
        if i =="c" or i == "d":
            lis.append(kwargs[i])

    def wrapper(*lis):
        print(*lis)

    wrapper(*lis)

function(**argument)
'''

argument = dict(
    b=2,
    c=0,
    d=4
)

def function(**kwargs):
    def wrapper(c, d):
        print(c, d)

    del kwargs["a"]
    del kwargs["b"]

    return wrapper(**kwargs)

wrap = function(a=0, b=1, c=2, d=3)
