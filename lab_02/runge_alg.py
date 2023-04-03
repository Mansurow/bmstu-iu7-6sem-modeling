import numpy as np
import matplotlib.pyplot as plt

from pandas.core.frame import DataFrame
from dataclasses import dataclass

@dataclass
class Lab:
    # Значения для рассчетов
    R: float    # Радиус трубки
    le: float   # Расстояние между электродами лампы
    Lk: float   # Индуктивность
    Ck: float   # Емкость конденсатора
    Rk: float   # Сопротивление
    Uco: float  # Напряжение на конденсаторе в начальный момент времени
    I0: float   # Сила тока в цепи в начальный момент времени t = 0
    Tw: float   # Температура

    # Введено из замечания
    pd: float
    Imax: float  # Максимальный ток

    # Значения для сходимости
    H = 1e-7  # Точность
    STEP = 1 / 30

    # Массивы данных
    iArr: list
    tkArr: list
    t0Arr: list
    sigmaArr: list
    mArr: list

    # Результат полученных методами Эйлера(Рунге-Кутта 1 порядка), метод Рунге-Кутта 2, 3 и 4 порядка
    iResArr = list()
    uResArr = list()
    # Значения полученных путем интегрирования (метод Трапеции)
    rResArr = list()
    zResArr = list()
    # Значения полученных путем интерполяции (Ньютона, Эрмита, сплайнами)
    tResArr = list()
    t0ResArr = list()
    mResArr = list()
    sigmaResArr = list()
    tkResArr = list()

    def __init__(self, i_t0_m: DataFrame, tk_sigma: DataFrame):
        self.R = 0.32  # см
        self.le = 12  # см
        self.Lk = 187 * (10**(-6))  # Гн
        self.Ck = 268 * (10**(-6))  # Ф
        self.Rk = 0.25  # Ом
        self.Uco = 1400  # B
        self.I0 = 0.3  # A
        self.Tw = 2000  # K

        self.t0Arr = i_t0_m.T0.array
        self.iArr = i_t0_m.I.array
        self.mArr = i_t0_m.m.array

        self.tkArr = tk_sigma.Tk.array
        self.sigmaArr = tk_sigma.sigma.array

    def dropResultArrs(self):
        self.tResArr.clear()
        self.iResArr.clear()
        self.uResArr.clear()
        self.rResArr.clear()
        self.t0ResArr.clear()
        self.sigmaResArr.clear()
        self.mResArr.clear()
        self.tkResArr.clear()

    # 1 уравнение из системы
    def dI_dt(self, i: float, u: float, r_p: float):
        return (u - (self.Rk + r_p) * i) / self.Lk

    # 2 уравнение из системы
    def dU_dt(self, i: float):
        return -(i / self.Ck)

    # функция нахождения T(z)
    def T(self, t0: float, z: float, m: float):
        return t0 + (self.Tw - t0) * np.power(z, m)

    # линейная интерполяция Ньютона
    def interpolate(self, arr1: list, arr2: list, value: float):
        n = len(arr1)
        j = 0

        if value < arr1[0]:
            return arr2[0]
        # elif value > arr1[n - 1]:
        #    return arr2[n - 1]

        while True:
            if arr1[j] > value or j == n - 2:
                break
            j += 1
        j -= 1

        if j < n - 1:
            dx = arr1[j + 1] - arr1[j]
            di = value - arr1[j]
            res = arr2[j] + ((arr2[j + 1] - arr2[j]) * di / dx)
            # print(j, i, i_arr[j+1], i_arr[j])
        else:
            dx = arr1[n - 1] - arr1[n - 2]
            di = value - arr1[n - 1]
            res = arr2[n - 2] + ((arr2[n - 1] - arr2[n - 2]) * di / dx)

        return res

    # интегрирование методом трапеции
    def rectangle_rule(self, func, a: float, b: float, h: float, i:float):
        s = (self.sigma_Т(a, i) + self.sigma_Т(b, i)) / 2
        while a < b + h:
            a += h
            s += self.sigma_Т(a, i)
        return s * h

    # f(z) = sigma(T(Z)) * z
    def sigma_Т(self, z: float, i: float):
        t0 = self.interpolate(self.iArr, self.t0Arr, i)
        m = self.interpolate(self.iArr, self.mArr, i)

        tk = self.T(t0, z, m)
        sigma = self.interpolate(self.tkArr, self.sigmaArr, tk)

        self.tkResArr.append(tk)
        self.zResArr.append(z)
        self.sigmaResArr.append(sigma)

        return sigma * z

    # нахождение Rp(I)
    def R_p(self, i):
        h = self.STEP
        z = 0
        z_max = 1

        r = self.le / (2 * np.pi * self.R * self.R *
                       self.rectangle_rule(self.sigma_Т, z, z_max, h, i))

        return r

    def euler(self, t=0, t_max=0.001, with_r=True):
        i_n = self.I0
        u_n = self.Uco
        t_n = t

        t0 = self.interpolate(self.iArr, self.t0Arr, i_n)
        m = self.interpolate(self.iArr, self.mArr, i_n)

        self.t0ResArr.append(t0)
        self.mResArr.append(m)
        self.tResArr.append(t_n)
        self.iResArr.append(i_n)
        self.uResArr.append(u_n)
        if with_r:
            r_n = self.R_p(i_n)
            self.rResArr.append(r_n)

        while t_n < t_max:
            if with_r:
                k1 = self.H * self.dI_dt(i_n, u_n, -self.Rk)
            else:
                k1 = self.H * self.dI_dt(i_n, u_n, -self.Rk)  # (Rk + Rp) == 0
            q1 = self.H * self.dU_dt(i_n)

            t_n = t_n + self.H
            i_n = i_n + k1
            u_n = u_n + q1
            if with_r:
                r_n = self.R_p(i_n)
                self.rResArr.append(r_n)

            t0 = self.interpolate(self.iArr, self.t0Arr, i_n)
            m = self.interpolate(self.iArr, self.mArr, i_n)

            self.t0ResArr.append(t0)
            self.mResArr.append(m)
            self.tResArr.append(t_n)
            self.iResArr.append(i_n)
            self.uResArr.append(u_n)

    def runge2(self, t=0, t_max=0.01, beta=1/2, with_r=True):
        i_n = self.I0
        u_n = self.Uco
        t_n = t
        t0 = self.interpolate(self.iArr, self.t0Arr, i_n)
        m = self.interpolate(self.iArr, self.mArr, i_n)

        self.t0ResArr.append(t0)
        self.mResArr.append(m)
        self.tResArr.append(t_n)
        self.iResArr.append(i_n)
        self.uResArr.append(u_n)
        if with_r:
            r_n = self.R_p(i_n)
            self.rResArr.append(r_n)

        while t_n < t_max:

            if with_r:
                k1 = self.H * self.dI_dt(i_n, u_n, r_n)
            else:
                k1 = self.H * self.dI_dt(i_n, u_n, -self.Rk)
            q1 = self.H * self.dU_dt(i_n)

            if with_r:
                k2 = self.H * self.dI_dt(i_n + k1 / (2 * beta), u_n + q1 / (2 * beta), self.R_p(i_n + k1 / (2 * beta)))
            else:
                k2 = self.H * self.dI_dt(i_n + k1 / (2 * beta), u_n + q1 / (2 * beta), -self.Rk)
            q2 = self.H * self.dU_dt(i_n + k1 / (2 * beta))

            t_n = t_n + self.H
            i_n = i_n + (1 - beta) * k1 + beta * k2
            u_n = u_n + (1 - beta) * q1 + beta * q2

            if with_r:
                r_n = self.R_p(i_n)
                self.rResArr.append(r_n)

            t0 = self.interpolate(self.iArr, self.t0Arr, i_n)
            m = self.interpolate(self.iArr, self.mArr, i_n)

            self.t0ResArr.append(t0)
            self.mResArr.append(m)
            self.tResArr.append(t_n)
            self.iResArr.append(i_n)
            self.uResArr.append(u_n)

    def runge4(self, t=0, t_max=0.01,  with_r=True):
        i_n = self.I0
        u_n = self.Uco
        t_n = t

        t0 = self.interpolate(self.iArr, self.t0Arr, i_n)
        m = self.interpolate(self.iArr, self.mArr, i_n)

        self.t0ResArr.append(t0)
        self.mResArr.append(m)
        self.tResArr.append(t_n)
        self.iResArr.append(i_n)
        self.uResArr.append(u_n)
        if with_r:
            r_k1 = self.R_p(i_n)
            self.rResArr.append(r_k1)

        while t_n < t_max:
            if with_r:
                k1 = self.H * self.dI_dt(i_n, u_n, r_k1)
                q1 = self.H * self.dU_dt(i_n)
                k2 = self.H * self.dI_dt(i_n + k1 / 2, u_n + q1 / 2, self.R_p(i_n + k1 / 2))
                q2 = self.H * self.dU_dt(i_n + k1 / 2)
                k3 = self.H * self.dI_dt(i_n + k2 / 2, u_n + q2 / 2, self.R_p(i_n + k2 / 2))
                q3 = self.H * self.dU_dt(i_n + k2 / 2)
                k4 = self.H * self.dI_dt(i_n + k3, u_n + q3, self.R_p(i_n + k3))

            else:
                k1 = self.H * self.dI_dt(i_n, u_n, -self.Rk)
                q1 = self.H * self.dU_dt(i_n)
                k2 = self.H * self.dI_dt(i_n + k1 / 2, u_n + q1 / 2, -self.Rk)
                q2 = self.H * self.dU_dt(i_n + k1 / 2)
                k3 = self.H * self.dI_dt(i_n + k2 / 2, u_n + q2 / 2, -self.Rk)
                q3 = self.H * self.dU_dt(i_n + k2 / 2)
                k4 = self.H * self.dI_dt(i_n + k3, u_n + q3, -self.Rk)
            q4 = self.H * self.dU_dt(i_n + k3)

            t_n = t_n + self.H
            i_n = i_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            u_n = u_n + (q1 + 2 * q2 + 2 * q3 + q4) / 6

            if with_r:
                r_k1 = self.R_p(i_n)
                self.rResArr.append(r_k1)

            t0 = self.interpolate(self.iArr, self.t0Arr, i_n)
            m = self.interpolate(self.iArr, self.mArr, i_n)

            self.t0ResArr.append(t0)
            self.mResArr.append(m)
            self.tResArr.append(t_n)
            self.iResArr.append(i_n)
            self.uResArr.append(u_n)

    def drawGraf(self, task=-1, text=""):
        print('Task' + str(task) + ": " + text)
        if task == 3:
            fig, axs = plt.subplots(1, 2, figsize=(10, 5))
            fig.suptitle('Task' + str(task) + ": " + text, fontsize=14, fontweight='bold')

            # Отрисовка графика функция I(t)
            axs[0].set_title("I(t)", fontsize=16, fontweight='bold')
            axs[0].set_xlabel('I')
            axs[0].set_ylabel('t')
            axs[0].plot(self.tResArr, self.iResArr)

            # Отрисовка графика функция U(t)
            axs[1].set_title("U(t)", fontsize=16, fontweight='bold')
            axs[1].set_xlabel('U')
            axs[1].set_ylabel('t')
            axs[1].plot(self.tResArr, self.uResArr)
        elif task == 2:
            fig, axs = plt.subplots(2, 3, figsize=(20, 15))
            fig.suptitle('Task' + str(task) + ": " + text, fontsize=14, fontweight='bold')

            # Отрисовка графика функция I(t)
            axs[0, 0].set_title("I(t)", fontsize=14, fontweight='bold')
            axs[0, 0].set_xlabel('t')
            axs[0, 0].set_ylabel('I')
            axs[0, 0].plot(self.tResArr, self.iResArr)

            # Отрисовка графика функция U(t)
            axs[0, 1].set_title("U(t)", fontsize=14, fontweight='bold')
            axs[0, 1].set_xlabel('t')
            axs[0, 1].set_ylabel('U')
            axs[0, 1].plot(self.tResArr, self.uResArr)

            # Отрисовка графика функция Rp(t)
            axs[0, 2].set_title("Rp(t)", fontsize=14, fontweight='bold')
            axs[0, 2].set_xlabel('t')
            axs[0, 2].set_ylabel('R')
            axs[0, 2].plot(self.tResArr, self.rResArr)

            # Отрисовка графика функция I(t) * Rp(t)
            i_mul_rResArr = [self.iResArr[i] * self.rResArr[i] for i in range(len(self.rResArr))]
            axs[1, 0].set_title("I(t) * Rp(t)", fontsize=14, fontweight='bold')
            axs[1, 0].set_xlabel('t')
            axs[1, 0].set_ylabel('I * Rp')
            axs[1, 0].plot(self.tResArr, i_mul_rResArr)

            # Отрисовка графика функция T0(t)
            axs[1, 1].set_title("T0(t)", fontsize=14, fontweight='bold')
            axs[1, 1].set_xlabel('t')
            axs[1, 1].set_ylabel('T0')
            axs[1, 1].plot(self.tResArr, self.t0ResArr)

            # Отрисовка графика функция sigma(t)
            axs[1, 2].set_title("sigma(t)", fontsize=14, fontweight='bold')
            axs[1, 2].set_xlabel('t')
            axs[1, 2].set_ylabel('sigma')
            self.tkResArr.sort()
            self.sigmaResArr.sort()
            axs[1, 2].plot(self.tkResArr, self.sigmaResArr)

        plt.show()

    def print_params(self):
        print("--------Params---------")
        print("1. R =", self.R)
        print("2. le =", self.le)
        print("3. Lk =", self.Lk)
        print("4. Ck =", self.Ck)
        print("5. Rk =", self.Rk)
        print("6. Uco =", self.Uco)
        print("7. I0 =", self.I0)
        print("8. Tw =", self.Tw)

        print("9. Длительность имульса =", self.pd)
        print("10. Imax =", self.Imax)
