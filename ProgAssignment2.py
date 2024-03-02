def neville(x, y, target_x):
    n = len(x)
    q = [[0] * n for _ in range(n)]

    for i in range(n):
        q[i][0] = y[i]

    for i in range(1, n):
        for j in range(1, i + 1):
            q[i][j] = ((target_x - x[i - j]) * q[i][j - 1] - (target_x - x[i]) * q[i - 1][j - 1]) / (x[i] - x[i - j])

    return q[n - 1][n - 1]

target_x = 3.7
x = [3.6, 3.8, 3.9]
y = [1.675, 1.436, 1.318]

res = neville(x, y, target_x)
print(f"Answer: {res}")

#Q2---------------------------------------------

def divided_difference(x, fx):
    if len(x) == 1:
        return fx[0]
        
    return (divided_difference(x[1:], fx[1:]) - divided_difference(x[:-1], fx[:-1])) / (x[-1] - x[0])

def newton_forward_coefficients(x, fx):
    n = len(x)
    coefficients = [fx[0]]
    
    for i in range(1, n):
        coefficients.append(divided_difference(x[:i+1], fx[:i+1]))
        
    return coefficients

def newton_forward_polynomial(x, fx, degree):
    coefficients = newton_forward_coefficients(x, fx)
    result = coefficients[0]
    term = 1
    degree += 1
    
    for i in range(1, degree):
        term *= (degree - i) / i * (7.4 - 7.2)
        result += term * coefficients[i]
        
    return result

#Data
x = [7.2, 7.4, 7.5, 7.6]
fx = [23.5492, 25.3913, 26.8224, 27.4589]

#Calculating degree 1, 2, and 3 approximations
for degree in range(1, 4):
    p = newton_forward_polynomial(x, fx, degree)
    
    print(f"Degree {degree}:", p)
    
 #Q3-----------------------------------------------

import math

def nfp2(x, fx, degree, table):
    h = x[1] - x[0]
    u = (7.3 - x[0]) / h
    result = fx[0]

    for i in range(1, degree + 1):
        term = table[0][i] / math.factorial(i)
        for j in range(i):
            term *= (u - j)
        result += term

    return result


#Data
x = [7.2, 7.4, 7.5, 7.6]
n = len(x)
fx = [23.5492, 25.3913, 26.8224, 27.4589]

#Calculating FDT
table = [[0 for _ in range(n)] for _ in range(n)]
for i in range(n):
    table[i][0] = fx[i]

for i in range(1, n):
    for j in range(n - i):
        table[j][i] = table[j + 1][i - 1] - table[j][i - 1]

#Approximating f(7.3) degree 1
degree = 1
fs = nfp2(x, fx, degree, table)
print(f"f(7.3)= ", fs)

##Q4-------------------------------

import numpy as np

#Data
x = np.array([3.6, 3.8, 3.9])
n = len(x) * 2
f_x = np.array([1.675, 1.436, 1.318])
f_prime_x = np.array([-1.195, -1.188, -1.182])

#Creating Matrix
matrix = np.zeros((n, n + 1))


for i in range(len(x)):
    matrix[i * 2, 0] = x[i]
    matrix[i * 2 + 1, 0] = x[i]
    matrix[i * 2, 1] = f_x[i]
    matrix[i * 2 + 1, 1] = f_x[i]
    matrix[i * 2, 2] = f_prime_x[i]
    matrix[i * 2 + 1, 2] = f_prime_x[i]


for i in range(3, n + 1):
    for j in range(n - i):
        if matrix[j + 1, 0] == matrix[j, 0]:
            matrix[j, i] = matrix[j + 1, i - 1]
        else:
            matrix[j, i] = (matrix[j + 1, i - 1] - matrix[j, i - 1]) / (matrix[j + 1, 0] - matrix[j, 0])


np.set_printoptions(precision = 8, suppress = True)
print("Hermite Matrix:")
print(matrix)

#Q5------------------------------------------

#Data
x = np.array([2, 5, 8, 10])
n = len(x)
f_x = np.array([3, 5, 7, 9])

#Making maxtix A and vector b
m = np.zeros((n, n))
b = np.zeros(n)


#Cubic Spline
m[0, 0] = 1
m[n - 1, n - 1] = 1

for i in range(1, n - 1):
    h_i = x[i] - x[i - 1]
    h_i_plus_1 = x[i + 1] - x[i]
    
    m[i, i - 1] = h_i
    m[i, i] = 2 * (h_i + h_i_plus_1)
    m[i, i + 1] = h_i_plus_1
    
    b[i] = 6 * ((f_x[i + 1] - f_x[i]) / h_i_plus_1 - (f_x[i] - f_x[i - 1]) / h_i)

#Solving for x
cx = np.linalg.solve(m, b)


a = f_x.copy()
b = np.zeros(n - 1)
d = np.zeros(n - 1)

for i in range(n - 1):
    h_i = x[i + 1] - x[i]
    b[i] = (a[i + 1] - a[i]) / h_i - h_i * (2 * cx[i] + cx[i + 1]) / 6
    d[i] = (cx[i + 1] - cx[i]) / (6 * h_i)


np.set_printoptions(precision=4, suppress=True)
print("Matrix A:")
print(m)
print("\nVector b:")
print(b)
print("\nVector x:")
print(cx)