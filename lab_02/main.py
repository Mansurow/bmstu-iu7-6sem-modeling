import pandas as pd
import matplotlib.pyplot as plt
from runge_alg import Lab

# Задана система электротехнических уравнений, описывающих разрядный контур,
# включающий постоянное активное сопротивление Rk, нелинейное сопротивление
# R (I) p, зависящее от тока I , индуктивность Lk и емкость Ck

PD_without_r = 0.01
PD = 600e-6  # Длительность импульса
Imax = 800

def print_menu():
    print("-------------------Menu-----------------------")
    print("1 - Euler")
    print("2 - Runge 2nd order")
    print("3 - Runge 4th order")
    print("4 - Euler Rk+Rp=0")
    print("5 - Runge 2nd Rk+Rp=0")
    print("6 - Runge 4nd Rk+Rp=0")
    print("7 - Runge 4nd Rk=200")
    print("0 - Exit")
    print("---------------------------------------------")

def main():
    table1 = pd.read_csv("./data/table1.csv")
    table2 = pd.read_csv("./data/table2.csv")

    print("-------------------Таблица 1-----------------")
    print(table2)
    print("---------------------------------------------")
    print("-------------------Таблица 2-----------------")
    print(table1)
    print("---------------------------------------------")

    f = Lab(table1, table2)
    choice = -1
    while choice != 0:
        print_menu()
        choice = int(input("Жду?:"))
        if choice == 1:
            f.dropResultArrs()
            f.euler(t=0, t_max=PD)
            f.drawGraf(task=2, text="Метод Эйлера")
        elif choice == 2:
            f.dropResultArrs()
            f.runge2(t=0, t_max=PD)
            f.drawGraf(task=2, text="Метод Рунге-Кутта 2 порядка")
        elif choice == 3:
            f.dropResultArrs()
            f.runge4(t=0, t_max=PD)
            f.drawGraf(task=2, text="Метод Рунге-Кутта 4 порядка")
        elif choice == 4:
            f.dropResultArrs()
            f.euler(t=0, t_max=PD_without_r, with_r=False)
            f.drawGraf(task=3, text="Метод Рунге-Кутта 4 порядка")
        elif choice == 5:
            f.dropResultArrs()
            f.runge2(t=0, t_max=PD_without_r, with_r=False)
            f.drawGraf(task=3, text="Метод Рунге-Кутта 4 порядка")
        elif choice == 6:
            f.dropResultArrs()
            f.runge4(t=0, t_max=PD_without_r, with_r=False)
            f.drawGraf(task=3, text="Метод Рунге-Кутта 4 порядка")
        elif choice == 7:
            f.Rk = 200
            print(f.Rk)
            f.runge4(0, 20e-6, with_r=True)

            plt.plot(f.tResArr, f.iResArr, "green", label="Runge 4")
            plt.legend(loc='best')
            plt.title("I(t)")
            plt.show()
            f.Rk = 0.25
if __name__ == '__main__':
   main()

