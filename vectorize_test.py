import numpy as np
def f(x1,x2):
    if x1 == 0:
        return x2 + 1
    else:
        return x1 + x2

x=np.array([1, 0, 3, 4])
y=np.array([2, 4, 5, 6])
#z=f(x, y)
z= np.vectorize(f)(x, y)
print(z)  # Output: [3 5 8 10]
# Note: The output for the vectorized function will be an array where each element is computed