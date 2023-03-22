import numpy as np
import matplotlib.pyplot as plt

a = 0  # left border
b = 1  # right border


def function(x):
    return x ** 4 - 2.9 * x ** 3 + 6.5 * x ** 2 - 7 * x - 5.4


value_of_integral = -7.258333333333335  # value of integrate function in [a,b] for checking


def method(left, right, section, epsilon):
    mas_error_1 = []
    mas_error_2 = []
    mas_section = []
    value = []
    summa_2 = 1e+5
    for k in range(15):
        sec = 2 ** k
        summa_1 = summa_2
        summa_2 = (function(left) + function(right)) / 2
        step = (abs(left - right)) / sec
        for j in range(1, sec):
            dx_1 = left + j * step
            summa_2 += function(dx_1)
        summa_2 = step * summa_2
        abs_error = (abs(summa_2 - summa_1)) / 3
        abs_error_2 = abs(summa_2 - value_of_integral)
        mas_error_1.insert(k, abs_error)
        mas_error_2.insert(k, abs_error_2)
        mas_section.insert(k, sec)
        value.insert(k, summa_2)
        if abs_error < epsilon or sec >= section:
            return mas_error_1, mas_error_2, abs_error, sec, mas_section
    return mas_error_1, mas_error_2, abs_error, sec, mas_section


# --------------------------------------------------graphs--------------------------------------------------------------
def extra_work():
    final_result = method(0, 1, 1024, 1e-7)
    res_1 = final_result[0]
    mas_sec = final_result[4]
    return mas_sec, res_1


def method_gauss(left, right, sec, epsilon):
    mas_error_1 = []
    mas_section = []
    mas_x = [-0.906180, -0.538469, 0.0, 0.538469, 0.906180]
    mas_a = [0.236927, 0.478629, 0.568889, 0.478629, 0.236927]
    sum_2 = 1e+5
    for i in range(15):
        sum_1 = sum_2
        node = 2 ** i + 1
        section = (abs(left - right)) / (node - 1)
        sum_2 = 0
        for j in range(1, node):
            for k in range(5):
                sum_2 += mas_a[k] * function(
                    (2 * left + (2 * j - 1) * section) / 2 + (section * mas_x[k]) / 2)
        sum_2 *= section / 2
        abs_error = (abs(sum_2 - sum_1)) / (2 ** (2 * 5) - 1)
        mas_error_1.insert(i, abs_error)
        mas_section.insert(i, node - 1)
        if abs_error < epsilon or (node - 1) >= sec:
            return mas_error_1, abs_error, node - 1, mas_section
    return mas_error_1, abs_error, node - 1, mas_section


final_result_1 = method_gauss(a, b, 1024, 1e-16)
mas_sec = final_result_1[3]
res_1_1 = final_result_1[0]

plt.figure(1)
plt.grid()
plt.yscale('log')
plt.xlabel('number of sections')
plt.ylabel('absolute error')
plt.title('Graph №1')
plt.plot(mas_sec[1:], res_1_1[1:], label="Runge criterion: 2 nodes")
plt.legend()


def graph_2():
    mas_tolerance = np.zeros((1, 7))[0]
    for j in range(7, 14):
        mas_tolerance[j - 7] = 10 ** (-j)
    mas_error = np.zeros((1, 7))[0]
    for j in range(7):
        mas_error[j] = method_gauss(a, b, 1024, mas_tolerance[j])[1]
    final = mas_tolerance, mas_error
    return final


plt.figure(2)
plt.grid()
plt.yscale('log')
plt.xscale('log')
plt.xlabel('specified accuracy')
plt.ylabel('experimental error')
plt.title('Graph №2')
plt.plot(graph_2()[0], graph_2()[0], label='best result')
plt.plot(graph_2()[0], graph_2()[1], label='real result: 2 nodes')
plt.legend()


def graph_3():
    mas_tolerance = np.zeros((1, 7))[0]
    for j in range(7, 14):
        mas_tolerance[j - 7] = 10 ** (-j)
    mas_n = np.zeros((1, 7))[0]
    for j in range(7):
        mas_n[j] = method_gauss(a, b, 1024, mas_tolerance[j])[2]
    final = mas_tolerance, mas_n
    return final


plt.figure(3)
plt.grid()
plt.xscale('log')
plt.xlabel('required accuracy')
plt.ylabel('required number of sections')
plt.title('Graph №3')
plt.plot(graph_3()[0], graph_3()[1], label='2 nodes')
plt.legend()


def graph_4():
    mas_error = []
    value = method_gauss(a, b, 1024, 1e-16)[0]
    section = len(value)
    for j in range(1, section):
        mas_error.insert(j, (value[1] / value[j]))
    return mas_error


plt.figure(4)
plt.grid()
plt.yscale('log')
plt.xlabel('number of sections')
plt.ylabel('absolute error')
plt.title('Graph №4')
x_5_1 = np.linspace(0, 1024, 10000)
y_5_1 = x_5_1 ** 7
plt.plot(mas_sec[1:], graph_4(), '--', label='2 nodes')
plt.plot(x_5_1, y_5_1, label='x^7')
plt.legend()

plt.figure(5)
plt.grid()
plt.yscale('log')
plt.xlabel('number of sections')
plt.ylabel('absolute error')
plt.title('Graph №5')
res_5_1 = extra_work()[0]
res_5_2 = extra_work()[1]
plt.plot(mas_sec[1:], res_1_1[1:], label="Runge criterion: method Gauss")
plt.plot(res_5_1[1:], res_5_2[1:], label="Runge criterion: method trapezoid")
plt.legend()
plt.show()
