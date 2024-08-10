import parameters as prm
from matplotlib import pyplot as plt


def generate_grid(start, stop, step):
    grid = [start]
    x = start
    while x < stop:
        x += step
        grid.append(x)
    if grid[-1] > stop:
        grid[-1] = stop
    return grid


def generate_table_value(table_fn, point):
    lt = len(table_fn)
    i = 0
    while i < lt - 2 and table_fn[i + 1][0] < point:
        i += 1
    dX = table_fn[i + 1][0] - table_fn[i][0]
    dY = table_fn[i + 1][1] - table_fn[i][1]
    k = dY / dX
    b = table_fn[i][1]
    dx = point - table_fn[i][0]
    dy = k * dx
    return b + dy


def generate_table_values(table_fn, points):
    return list(map(lambda x: generate_table_value(table_fn, x), points))


def generate_func_values(func, points):
    return list(map(lambda x: func(x), points))


def generate_abcd(r_points, T_values, prev_T_values, time_delta, time_value=0,
                  F_const_mode=0, F_frequency=0):
    h = r_points[1] - r_points[0]
    N = len(r_points)
    lambdas = generate_table_values(prm.LAMBDA, T_values)
    ks = generate_table_values(prm.K, T_values)
    if F_const_mode:
        F_ = prm.F_const(time_value)
    elif F_frequency > 0:
        F_length = 1 / F_frequency * 10 ** 6  # с в мкс
        transform_time_value = time_value % F_length
        F_ = prm.F(transform_time_value)
    else:
        F_ = prm.F(time_value)
    res = [
        [
            0,
            -lambdas[0] / h,
            -lambdas[0] / h,
            -F_
        ]
    ]
    for n in range(1, N - 1):
        a = lambdas[n] / (prm.R * r_points[n]) + (
                lambdas[n + 1] - lambdas[n]) / h
        b = lambdas[n]
        if F_const_mode:
            c = prm.c_NEAR_ZERO
        else:
            c = prm.a2 + prm.b2 * T_values[n] ** prm.m2 - prm.c2 / T_values[
                n] ** 2
        f = 4 * ks[n] * prm.np ** 2 * prm.sigma * (
                T_values[n] ** 4 - prm.T0 ** 4
        )
        K1 = b / h ** 2
        K3 = K1 + a / h
        K2 = c / (0.5 * time_delta) + K1 + K3
        K4 = c * prev_T_values[n] / (0.5 * time_delta) - f
        res.append([K1, K2, K3, K4])
    res.append([
        lambdas[-1] / h,
        lambdas[-1] / h,
        0,
        -prm.ALPHA * (T_values[-1] - prm.T0)
    ])
    return res


def do_sweep(abcd, points):
    N = len(points)
    ksi = [0 for _ in range(N)]
    eta = [0 for _ in range(N)]

    for i in range(1, N):
        den = (abcd[i - 1][1] - abcd[i - 1][0] * ksi[i - 1])
        ksi[i] = abcd[i - 1][2] / den
        eta[i] = (abcd[i - 1][3] + abcd[i - 1][0] * eta[i - 1]) / den

    y = [0 for _ in range(N)]
    y[N - 1] = (abcd[-1][3] + abcd[-1][0] * eta[-1]) / (
            abcd[-1][1] - abcd[-1][0] * ksi[-1])
    for i in range(N - 2, -1, -1):
        y[i] = y[i + 1] * ksi[i + 1] + eta[i + 1]
    return y


def solution(f_mode=0):
    points = generate_grid(prm.r0 / prm.R, 1, prm.h)
    N = len(points)
    prev_ts = [prm.T0 for _ in points]
    t = 0
    while t < prm.max_t:
        it = 0
        ts = prev_ts
        while it < prm.max_iters:
            it += 1
            abcd = generate_abcd(points, ts, prev_ts, prm.tau, t, f_mode)
            y = do_sweep(abcd, points)
            max_delta = abs(y[0] - ts[0]) / abs(ts[0])
            for i in range(1, N):
                max_delta = max(max_delta, abs(y[i] - ts[i]) / abs(ts[i]))
            if max_delta < prm.temp_eps:
                break
            ts = [(ts[i] + y[i]) / 2 for i in range(N)]
        yield t, points, ts, it
        t += prm.tau
        prev_ts = ts


def solution_impulse(f_mode=0, f_freq=0):
    points = generate_grid(prm.r0 / prm.R, 1, prm.h)
    N = len(points)
    prev_ts = [prm.T0 for _ in points]
    t = 0
    n_it = 0
    max_n_it = 10
    while t <= prm.max_t:  # /10:
        it = 0
        ts = prev_ts
        while it < prm.max_iters:
            it += 1
            abcd = generate_abcd(points, ts, prev_ts, prm.tau, t, f_mode,
                                 f_freq)
            y = do_sweep(abcd, points)
            max_delta = abs(y[0] - ts[0]) / abs(ts[0])
            for i in range(1, N):
                max_delta = max(max_delta, abs(y[i] - ts[i]) / abs(ts[i]))
            if max_delta < prm.temp_eps:
                break
            ts = [(ts[i] + y[i]) / 2 for i in range(N)]
        yield t, points, ts, it
        t += prm.tau
        prev_ts = ts


def T_z_t(F_const_mode=0):
    plt.figure(figsize=(6, 6))
    plt.xlabel("r0/R")
    plt.ylabel("T, K")
    if F_const_mode:
        plt.title(f"Результат работы программы для параметров\n"
                  f"F0={prm.F0},alpha={prm.ALPHA},h={prm.h},tau={prm.tau}.\n")
        max_it = 11
    else:
        plt.title(f"Результат работы программы для параметров\n"
                  f"Fmax={prm.F_MAX}, tmax={prm.t_MAX}, alpha={prm.ALPHA},h={prm.h},tau={prm.tau}.\n")
        max_it = prm.MAX_it
    plt.xlim(0.7, 1)
    # plt.ylim(0, 3000)
    it = 0
    for [_t, _points, _ts, _it] in solution(F_const_mode):
        it += 1
        if it % 10 == 1 and it < max_it:
            plt.plot(_points, _ts, label=f"t={_t:.6f},it={_it}")
    plt.legend()
    plt.show()


def T_0_t():
    plt.figure(figsize=(6, 6))
    plt.xlabel("t, мкс")
    plt.ylabel("T, K")
    plt.title(f"Функция T(0,t)\n"
              f"Fmax={prm.F_MAX}, tmax={prm.t_MAX}, alpha={prm.ALPHA},h={prm.h},tau={prm.tau}.\n")
    # plt.xlim(0.7, 1)
    # plt.ylim(0, 3000)
    it = 0
    x_arr = []
    y_arr = []
    for [_t, _points, _ts, _it] in solution():
        print(_t)
        x_arr.append(_t)
        y_arr.append(_ts[0])
    plt.plot(x_arr, y_arr)
    plt.show()


def T_0_t_impulse():
    plt.figure(figsize=(6, 6))
    plt.xlabel("t, мкс")
    plt.ylabel("T, K")
    plt.title(f"Функция T(0,t)\n"
              f"Fmax={prm.F_MAX}, tmax={prm.t_MAX}, alpha={prm.ALPHA},h={prm.h},tau={prm.tau},freq={prm.F_freq}Hz.\n")
    # plt.xlim(0.7, 1)
    # plt.ylim(0, 3000)
    it = 0
    x_arr = []
    y_arr = []
    for [_t, _points, _ts, _it] in solution_impulse(0, prm.F_freq):
        print(_t)
        x_arr.append(_t)
        y_arr.append(_ts[0])
    plt.plot(x_arr, y_arr)
    plt.show()


def main():
    #T_z_t(1)
    #T_0_t()
    T_0_t_impulse()


if __name__ == '__main__':
    main()
