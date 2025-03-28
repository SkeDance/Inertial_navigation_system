import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Константы
PI = math.pi
ECC = 0.081819190842622
g = 9.81523
U = 15 * (PI / 180) / 3600
R_Earth = 6.371e6
dt = 0.005

# Инициализация переменных
alignment_flag = False
t = 1
fi_0 = lambda_0 = fi = lambda_ = 0.0
ROLL_0 = PITCH_0 = YAW_0 = 0.0
current_altitude = 0.0

matrix_LL = np.identity(3)
Acc_matrix_BL = np.zeros((3, 1))
Gyro_matrix_BL = np.zeros((3, 1))
Acc_matrix_ENUp = np.zeros((3, 1))

VE = VN = 0.0
prev_aX = curr_aX = prev_aY = curr_aY = 0.0
integral_VN = integral_VE = 0.0

# Данные для графиков
time_points, pitch_points, roll_points, yaw_points = [], [], [], []
file_time, file_pitch, file_roll, file_yaw = [], [], [], []
file_fi, file_lambda, file_altitude = [], [], []
calc_fi, calc_lambda = [], []
deltas_roll, deltas_pitch, deltas_yaw = [], [], []
deltas_fi, deltas_lambda = [], []


def normalize_angle(degrees):
    degrees = math.fmod(degrees, 360)
    return degrees + 360 if degrees < 0 else degrees


def degrees_to_rads(degrees):
    return degrees * PI / 180


def rads_to_degrees(rads):
    return rads * 180 / PI


def body_to_local(A, B):
    return np.dot(A, B)


def update_orientation_matrix(LL, W_B, W_LL):
    LL_dt = np.dot(LL, W_B) - np.dot(W_LL, LL)
    return LL + LL_dt * dt


def normalize_matrix(C, step):
    if step % 2 == 0:
        return C / np.linalg.norm(C, axis=1)[:, np.newaxis]
    return C / np.linalg.norm(C, axis=0)


def calculate_orientation(matrix):
    C_0 = math.hypot(matrix[2][0], matrix[2][2])
    PITCH = round(rads_to_degrees(math.atan(matrix[2][1] / C_0)), 5)
    ROLL = round(-rads_to_degrees(math.atan2(matrix[2][0], matrix[2][2])), 5)
    YAW = round(normalize_angle(rads_to_degrees(math.atan2(matrix[0][1], matrix[1][1]))), 5)
    return PITCH, ROLL, YAW


def process_line(tokens):
    global alignment_flag, matrix_LL, t, fi, lambda_, current_altitude
    global VE, VN, prev_aX, curr_aX, prev_aY, curr_aY, integral_VN, integral_VE
    global fi_0, lambda_0, ROLL_0, PITCH_0, YAW_0
    global deltas_roll, deltas_pitch, deltas_yaw, deltas_fi, deltas_lambda

    try:
        time = round(float(tokens[0]), 5)
        Gyro_matrix_BL[:, 0] = [degrees_to_rads(float(x)) for x in tokens[1:4]]
        Acc_matrix_BL[:, 0] = [float(x) for x in tokens[4:7]]

        if len(tokens) >= 10:
            file_time.append(time)
            file_roll.append(round(float(tokens[7]), 5))
            file_pitch.append(round(float(tokens[8]), 5))
            file_yaw.append(round(float(tokens[9]), 5))

            if len(tokens) >= 15:
                file_fi.append(degrees_to_rads(float(tokens[12])))
                file_lambda.append(degrees_to_rads(float(tokens[13])))
                current_altitude = float(tokens[14])
                file_altitude.append(current_altitude)

        if not alignment_flag:
            if len(tokens) >= 15:
                ROLL_0 = round(float(tokens[7]), 5)
                PITCH_0 = round(float(tokens[8]), 5)
                YAW_0 = round(float(tokens[9]), 5)
                fi_0 = degrees_to_rads(float(tokens[12]))
                lambda_0 = degrees_to_rads(float(tokens[13]))

            cp = math.cos(degrees_to_rads(PITCH_0))
            cr = math.cos(degrees_to_rads(ROLL_0))
            sp = math.sin(degrees_to_rads(PITCH_0))
            sr = math.sin(degrees_to_rads(ROLL_0))
            cy = math.cos(degrees_to_rads(YAW_0))
            sy = math.sin(degrees_to_rads(YAW_0))

            matrix_LL = np.array([
                [cr * cy + sp * sr * sy, -cr * sy + sp * sr * cy, cp * sr],
                [cp * sy, cp * cy, sp],
                [sr * cy - sp * cr * sy, -sr * sy - sp * cr * cy, cp * cr]
            ])
            alignment_flag = True
            fi, lambda_ = fi_0, lambda_0

        Acc_matrix_ENUp = body_to_local(matrix_LL, Acc_matrix_BL)

        prev_aX, curr_aX = curr_aX, Acc_matrix_ENUp[0][0]
        VE = round(VE + (prev_aX + curr_aX) / 2 * dt, 8)

        prev_aY, curr_aY = curr_aY, Acc_matrix_ENUp[1][0]
        VN = round(VN + (prev_aY + curr_aY) / 2 * dt, 8)

        sin_fi = math.sin(fi)
        R_fi = round((R_Earth * (1 - ECC ** 2)) / (1 - ECC ** 2 * sin_fi ** 2) ** 1.5, 8)
        R_lambda = round(R_Earth / math.sqrt(1 - ECC ** 2 * sin_fi ** 2), 8)

        integral_VN = round(integral_VN + (VN / (R_fi + current_altitude)) * dt, 8)
        integral_VE = round(integral_VE + (VE / ((R_lambda + current_altitude) * math.cos(fi))) * dt, 8)
        fi = round(fi_0 + integral_VN, 8)
        lambda_ = round(lambda_0 + integral_VE, 8)

        calc_fi.append(fi)
        calc_lambda.append(lambda_)

        wE = round(-VN / (R_fi + current_altitude), 8)
        wN = round(VE / (R_lambda + current_altitude) + U * math.cos(fi), 8)
        wUp = round(VE * math.tan(fi) / (R_lambda + current_altitude) + U * math.sin(fi), 8)

        W_LL = np.array([
            [0, -wUp, wN],
            [wUp, 0, -wE],
            [-wN, wE, 0]
        ])

        W_B = np.array([
            [0, -Gyro_matrix_BL[2][0], Gyro_matrix_BL[1][0]],
            [Gyro_matrix_BL[2][0], 0, -Gyro_matrix_BL[0][0]],
            [-Gyro_matrix_BL[1][0], Gyro_matrix_BL[0][0], 0]
        ])

        matrix_LL = np.round(update_orientation_matrix(matrix_LL, W_B, W_LL), 8)
        matrix_LL = normalize_matrix(matrix_LL, t)
        t += 1

        PITCH, ROLL, YAW = calculate_orientation(matrix_LL)
        time_points.append(time)
        pitch_points.append(PITCH)
        roll_points.append(ROLL)
        yaw_points.append(YAW)

        if len(tokens) >= 10:
            delta_r = round(ROLL - file_roll[-1], 5)
            delta_p = round(PITCH - file_pitch[-1], 5)
            delta_y = round(YAW - file_yaw[-1], 5)

            if delta_y > 180:
                delta_y -= 360
            elif delta_y < -180:
                delta_y += 360

            deltas_roll.append(abs(delta_r))
            deltas_pitch.append(abs(delta_p))
            deltas_yaw.append(abs(delta_y))

            if len(file_fi) == len(calc_fi):
                delta_fi = rads_to_degrees(abs(calc_fi[-1] - file_fi[-1]))
                delta_lambda = rads_to_degrees(abs(calc_lambda[-1] - file_lambda[-1]))
                deltas_fi.append(delta_fi)
                deltas_lambda.append(delta_lambda)

            print(f"[t={time:.1f}s] Высота: {current_altitude:.1f}m | "
                  f"Roll Δ: {delta_r:+07.3f}° | Pitch Δ: {delta_p:+07.3f}° | Yaw Δ: {delta_y:+07.3f}°")

    except Exception as e:
        print(f"Ошибка обработки строки: {str(e)}")


def plot_orientation_results():
    plt.figure(figsize=(12, 10))
    plt.rcParams['axes.formatter.useoffset'] = False
    plt.rcParams['axes.formatter.limits'] = [-5, 5]

    for i, (title, calc, ref, ylabel) in enumerate(zip(
            ['Тангаж (Pitch)', 'Крен (Roll)', 'Курс (Yaw)'],
            [pitch_points, roll_points, yaw_points],
            [file_pitch, file_roll, file_yaw],
            ['Угол [°]', 'Угол [°]', 'Угол [°]']
    )):
        plt.subplot(3, 1, i + 1)
        plt.plot(time_points, calc, 'b-', linewidth=1.5, label='Расчёт')
        if ref and len(ref) == len(time_points):
            plt.plot(file_time, ref, 'r--', linewidth=1, label='Файл')
        plt.title(title)
        plt.ylabel(ylabel)
        plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
        plt.grid(True)
        plt.legend()

    plt.tight_layout()
    plt.savefig('orientation_comparison.png', dpi=300)
    plt.close()


def plot_coordinates_comparison():
    if not file_fi or not calc_fi:
        print("\nНедостаточно данных для сравнения координат")
        return

    plt.figure(figsize=(12, 8))

    # График широты
    plt.subplot(2, 1, 1)
    plt.plot(time_points, [rads_to_degrees(f) for f in calc_fi],
             'b-', label='Расчёт', linewidth=1.5)
    plt.plot(file_time[:len(file_fi)], [rads_to_degrees(f) for f in file_fi],
             'r--', label='Файл', linewidth=1)
    plt.title('Сравнение широты')
    plt.ylabel('Широта [°]')
    plt.grid(True)
    plt.legend()

    # График долготы
    plt.subplot(2, 1, 2)
    plt.plot(time_points, [rads_to_degrees(l) for l in calc_lambda],
             'b-', label='Расчёт', linewidth=1.5)
    plt.plot(file_time[:len(file_lambda)], [rads_to_degrees(l) for l in file_lambda],
             'r--', label='Файл', linewidth=1)
    plt.title('Сравнение долготы')
    plt.xlabel('Время [с]')
    plt.ylabel('Долгота [°]')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.savefig('coordinates_comparison.png', dpi=300)
    plt.close()


def plot_altitude():
    if not file_altitude:
        print("\nНет данных о высоте")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(file_time[:len(file_altitude)], file_altitude, 'g-', linewidth=2)
    plt.title('Высота от времени')
    plt.xlabel('Время [с]')
    plt.ylabel('Высота [м]')
    plt.grid(True)
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    plt.tight_layout()
    plt.savefig('altitude.png', dpi=300)
    plt.close()


def plot_coordinates_scatter():
    if not file_fi or not file_lambda:
        print("\nНедостаточно данных для построения графика координат")
        return

    plt.figure(figsize=(10, 8))

    file_lat = [rads_to_degrees(f) for f in file_fi]
    file_lon = [rads_to_degrees(l) for l in file_lambda]
    calc_lat = [rads_to_degrees(f) for f in calc_fi]
    calc_lon = [rads_to_degrees(l) for l in calc_lambda]

    plt.scatter(file_lon, file_lat, c='blue', s=10, label='Файл', alpha=0.7)
    plt.scatter(calc_lon, calc_lat, c='red', s=10, label='Расчёт', alpha=0.7)

    plt.title('Траектория движения')
    plt.xlabel('Долгота [°]')
    plt.ylabel('Широта [°]')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.gca().xaxis.set_major_formatter(plt.FormatStrFormatter('%.6f'))
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.6f'))
    plt.tight_layout()
    plt.savefig('coordinates_scatter.png', dpi=300)
    plt.close()


def print_max_deltas():
    print("\nМаксимальные расхождения:")
    if deltas_roll:
        print(f"Крен:    {max(deltas_roll):.5f}°")
    if deltas_pitch:
        print(f"Тангаж:  {max(deltas_pitch):.5f}°")
    if deltas_yaw:
        print(f"Курс:    {max(deltas_yaw):.5f}°")
    if deltas_fi:
        print(f"Широта:  {max(deltas_fi):.8f}°")
    if deltas_lambda:
        print(f"Долгота: {max(deltas_lambda):.8f}°")


def main():
    try:
        print("Запуск обработки данных...")
        with open("20170922_Cessna172_200Hz_Ref.txt", "r") as f:
            for line in f:
                tokens = line.strip().split()
                if len(tokens) >= 10:
                    process_line(tokens)

        print("\nГенерация графиков...")
        plot_orientation_results()
        plot_coordinates_comparison()
        plot_coordinates_scatter()
        plot_altitude()

        print_max_deltas()
        print("\nОбработка завершена. Графики сохранены в файлы:")
        print("- orientation_comparison.png")
        print("- coordinates_comparison.png")
        print("- coordinates_scatter.png")
        print("- altitude.png")

    except FileNotFoundError:
        print("Ошибка: файл с данными не найден!")
    except Exception as e:
        print(f"Критическая ошибка: {str(e)}")


if __name__ == "__main__":
    main()