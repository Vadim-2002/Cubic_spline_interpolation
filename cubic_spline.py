from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import math


def cubic_spline(x_i, y, x, h, m, quantity, i):
    S_3 = (((x_i[i] - x)**3 - (h**2)*(x_i[i] - x))*m[i - 1])/(6*h) + (((x - x_i[i - 1])**3 - (h**2)*(x - x_i[i - 1]))*m[i])/(6*h) + (((x_i[i] - x)*y[i - 1])/h) + ((((x - x_i[i - 1])*y[i]))/h)
    
    return S_3


def give_points(left_border, right_border, h, quantity):
    quantity += 1 #кол-во точек
    
    x = [i for i in range(quantity)]
    for i in range(0, quantity):
        x[i] = left_border + i*h

    y = [i for i in range(quantity)]
    for i in range(0, quantity):
        y[i] = np.sin(x[i])
    y[quantity - 1] = 0.0
    
    return x, y


def run_through_method(y, h, quantity):
    size = quantity + 1
    n = quantity
    
    a = [0 * size for i in range(size)]
    b = [0 * size for i in range(size)]
    c = [0 * size for i in range(size)]
    d = [0 * size for i in range(size)]
    m = [0 * size for i in range(size)]
    
    a[0] = 1
    b[0] = 0
    d[0] = 0
    c[n] = 0
    a[n] = 1
    d[n] = 0
    
    for i in range(1, n):
        a[i] = (h + h)/3
        c[i] = h/6
        b[i] = h/6
        d[i] = (y[i + 1] + y[i - 1] - 2*y[i]) / h

    abc = [[0] * size for i in range(size)]
    
    for i in range(0, size):
        for j in range(0, size):
            if(i == j):
                abc[i][j] = a[i]
            elif(i-1 == j):
                abc[i][j] = c[i]
            elif(i == j-1):
                abc[i][j] = b[i]
    
    lam = [0 * size for i in range(size)]
    u = [0 * size for i in range(size)]

    lam[0] = -b[0]/a[0]
    u[0] = d[0]/a[0]
    
    for i in range(1, size):
        lam[i] = -b[i]/(a[i] + c[i]*lam[i - 1])
        
    for i in range(1, size):
        u[i] = (d[i] - c[i]*u[i - 1])/(a[i] + c[i]*lam[i - 1])
    
    #решение системы:
    m[n] = u[n]
    for i in range(n - 1, -1, -1):
        m[i] = lam[i]*m[i + 1] + u[i]
    
    return m


def give_delta_max(left_border, right_border, quantity):
    h = (right_border - left_border) / quantity
    x = give_points(left_border, right_border, h, quantity)[0]
    y = give_points(left_border, right_border, h, quantity)[1]
    
    m = run_through_method(y, h, quantity)

    S = [0 * i for i in range(quantity + 1)]
    x_average = [0 * i for i in range(quantity + 1)]
    
    delta_max = -1
    for i in range(1, quantity + 1):
        x_average[i] = (x[i - 1] + x[i]) / 2
        
        S[i] = cubic_spline(x, y, x_average[i], h, m, quantity, i)
        
        dif = abs(S[i] - np.sin((x_average[i])))
        if (delta_max < dif):
            delta_max = dif
    
    return delta_max


def error_estimation(left_border, right_border):
    quantity = [0 * i for i in range(12)]    #кол-во отрезков
    delta_max = [0 * i for i in range(12)]   #максимальная погрешность для каждого отрезка
    delta_k = [0 * i for i in range(12)]     #отношение погрешности предыдущей строки к данной
    estimation = [0 * i for i in range(12)]  #оценка погрешности
    
    quantity[0] = 5
    for i in range(1, 12):
        quantity[i] = quantity[i - 1] * 2

    for i in range(0, 12):
        delta_max[i] = give_delta_max(left_border, right_border, quantity[i])

    for i in range(1, 12):
        delta_k[i] = delta_max[i - 1] / delta_max[i]
        estimation[i] = delta_max[i - 1] / (2**(round(math.log(delta_k[i], 2))))

    percentile_list = pd.DataFrame(
    {'n': quantity,
     'del_max': delta_max,
     'del_оц': estimation,
     'K': delta_k
    })
    print(percentile_list)


def draw_graph(x, y):
    plt.plot(x, y)   
    plt.ylabel("Y")   
    plt.xlabel("X")   
    plt.show()


if __name__ == '__main__':
    left_border = 0
    right_border = np.pi
    
    error_estimation(left_border, right_border)

    quantity = 50
    h = (right_border - left_border) / quantity
    x = give_points(left_border, right_border, h, quantity)[0]
    y = give_points(left_border, right_border, h, quantity)[1]
    
    m = run_through_method(y, h, quantity)

    S = [0 * i for i in range(quantity + 1)]
    x_average = [0 * i for i in range(quantity + 1)]
    
    for i in range(1, quantity + 1):
        x_average[i] = (x[i - 1] + x[i]) / 2
        
        S[i] = cubic_spline(x, y, x_average[i], h, m, quantity, i)
    
    draw_graph(x_average, S)
