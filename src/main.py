from thermo import cp, enthalpy, mixture_cp
from mixture import SPECIES_NAMES

def main():
    print("Расчёт свойств газов (NASA-коэффициенты)\n")

    for name in SPECIES_NAMES:
        print(f"{name:5s}  Cp(300K) = {cp(name, 300):8.2f}  Cp(1500K) = {cp(name, 1500):8.2f}")


    air = {"N2": 0.79, "O2": 0.21}
    print(f"\nВоздух Cp(500K) = {mixture_cp(air, 500):.2f} Дж/(кг·К)")


if __name__ == "__main__":
    main()
