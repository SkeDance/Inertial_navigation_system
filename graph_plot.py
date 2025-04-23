import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from math import radians, sin, cos, sqrt, atan2


def haversine(lat1, lon1, lat2, lon2):
    """
    Рассчитывает расстояние между двумя точками на сфере
    (формула гаверсинусов).
    Возвращает расстояние в километрах.
    """
    # Конвертация градусов в радианы
    lat1, lon1, lat2, lon2 = map(radians, [lat1, lon1, lat2, lon2])

    # Разницы координат
    dlat = lat2 - lat1
    dlon = lon2 - lon1

    # Формула гаверсинусов
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    # Средний радиус Земли в км
    R = 6371.0
    return R * c

# =============================================
# Загрузка и обработка данных
# =============================================

# Расчетные данные
calc_data = np.loadtxt('errors.txt', skiprows=1, delimiter='\t')
time_calc = calc_data[:, 0] - calc_data[0, 0]
lat_calc = calc_data[:, 4]
lon_calc = calc_data[:, 5]
roll_calc = calc_data[:, 2]
pitch_calc = calc_data[:, 1]
yaw_calc = calc_data[:, 3]
calc_val7 = calc_data[:, 6]  # 7-й столбец расчетных данных
calc_val8 = calc_data[:, 7]  # 8-й столбец расчетных данных

# Референсные данные
ref_data = np.loadtxt('20170922_Cessna172_200Hz_Ref.txt', delimiter=None)
time_ref = ref_data[:, 0] - ref_data[0, 0]
lat_ref = ref_data[:, 12]
lon_ref = ref_data[:, 13]
alt_ref = ref_data[:, 14]
roll_ref = ref_data[:, 7]
pitch_ref = ref_data[:, 8]
yaw_ref = ref_data[:, 9]
ref_val11 = ref_data[:, 10]  # 11-й столбец референсных данных
ref_val12 = ref_data[:, 11]  # 12-й столбец референсных данных

# =============================================
# Вычисление отклонений
# =============================================

# Интерполяция расчетных данных на временную сетку референсных данных
roll_calc_interp = np.interp(time_ref, time_calc, roll_calc)
pitch_calc_interp = np.interp(time_ref, time_calc, pitch_calc)
yaw_calc_interp = np.interp(time_ref, time_calc, yaw_calc)

# Вычисление разниц
diff_roll = roll_ref - roll_calc_interp
diff_pitch = pitch_ref - pitch_calc_interp
diff_yaw = yaw_ref - yaw_calc_interp
diff_ve = ref_val11 - calc_val7
diff_vn = ref_val12 - calc_val8

## Нахождение максимальных абсолютных отклонений
max_roll = np.max(np.abs(diff_roll))
max_pitch = np.max(np.abs(diff_pitch))
max_yaw = np.max(np.abs(diff_yaw)) - 359.95712
max_ve = np.max(np.abs(diff_ve))
max_vn = np.max(np.abs(diff_vn))

# Ошибки по курсу на первом и последнем тактах
yaw_error_start = diff_yaw[0]
yaw_error_end = diff_yaw[-1]

# Интерполяция расчетных координат
lat_calc_interp = np.interp(time_ref, time_calc, lat_calc)
lon_calc_interp = np.interp(time_ref, time_calc, lon_calc)

# Находим индекс 15-й минуты (900 секунд)
idx_15min = np.argmin(np.abs(time_ref - 900))

# Вычисляем ошибки координат
error_lat = lat_ref[idx_15min] - lat_calc_interp[idx_15min]
error_lon = lon_ref[idx_15min] - lon_calc_interp[idx_15min]

# Дополнительно: ошибка в метрах (1° ≈ 111 км)
lat_error_m = error_lat * 111000  # 1° широты ≈ 111 км
lon_error_m = error_lon * 111000 * np.cos(np.radians(lat_ref[idx_15min]))

# Находим индекс 15-й минуты (900 секунд)
idx_15min_ref = np.argmin(np.abs(time_ref - 900))
idx_15min_calc = np.argmin(np.abs(time_calc - 900))

# Координаты точек
lat_ref_pt = lat_ref[idx_15min_ref]
lon_ref_pt = lon_ref[idx_15min_ref]
lat_calc_pt = lat_calc[idx_15min_calc]
lon_calc_pt = lon_calc[idx_15min_calc]

# Расчет расстояния
distance_km = haversine(lat_ref_pt, lon_ref_pt, lat_calc_pt, lon_calc_pt)

# Вывод результатов
print("\nМаксимальные отклонения:")
print(f"Крен: {max_roll:.4f} град")
print(f"Тангаж: {max_pitch:.4f} град")
print(f"Курс_max: {max_yaw:.4f} град")
print(f"Курс_0: {yaw_error_start:.4f} град")
print(f"Курс_end: {yaw_error_end:.4f} град\n")
print("\nОшибка по координатам на 15-й минуте (900 сек):")
print(f"Широта: {error_lat:.8f} град | {lat_error_m:.2f} м")
print(f"Долгота: {error_lon:.8f} град | {lon_error_m:.2f} м")
print("\nРасстояние между точками на 15-й минуте:")
print(f"Референсная точка:  {lat_ref_pt:.6f}°, {lon_ref_pt:.6f}°")
print(f"Расчетная точка:    {lat_calc_pt:.6f}°, {lon_calc_pt:.6f}°")
print(f"Расстояние: {distance_km:.3f} км")
print(f"VE:  {max_ve: .4f} м/с")
print(f"VN: {max_vn: .4f} м/с\n")



# =============================================
# График траекторий в координатах
# =============================================
plt.figure(figsize=(14, 8))
plt.plot(lon_ref, lat_ref, 'r-', label='Истинная траектория', linewidth=1.5, alpha=0.8)
plt.plot(lon_calc, lat_calc, 'b--', label='Рассчитанная траектория', linewidth=1.5, alpha=0.8)

# Параметры зоны безопасности
radius_nm = 2
radius_km = radius_nm * 1.852

# Референсные данные: точка + эллипс
ref_15min_idx = np.argmin(np.abs(time_ref - 900))
lat_ref_pt = lat_ref[ref_15min_idx]
lon_ref_pt = lon_ref[ref_15min_idx]

delta_lat_ref = radius_km / 111.0
delta_lon_ref = radius_km / (111.0 * np.cos(np.radians(lat_ref_pt)))

plt.plot(lon_ref_pt, lat_ref_pt, 'o',
         markersize=12,
         markerfacecolor='none',
         markeredgecolor='magenta',
         markeredgewidth=2,
         label=f'Истинная траектория: 15 мин')

ellipse_ref = Ellipse(
    (lon_ref_pt, lat_ref_pt),
    width=2*delta_lon_ref,
    height=2*delta_lat_ref,
    edgecolor='magenta',
    linestyle='--',
    facecolor='none',
    linewidth=1.5,
    alpha=0.7
)
plt.gca().add_patch(ellipse_ref)

# Расчетные данные: точка
calc_15min_idx = np.argmin(np.abs(time_calc - 900))
lat_calc_pt = lat_calc[calc_15min_idx]
lon_calc_pt = lon_calc[calc_15min_idx]

plt.plot(lon_calc_pt, lat_calc_pt, 'o',
         markersize=5,
         markerfacecolor='none',
         markeredgecolor='lime',
         markeredgewidth=2,
         label='Расчёт: 15 мин')

plt.xlabel('Долгота [град]', fontsize=12)
plt.ylabel('Широта [град]', fontsize=12)
plt.title('Сравнение траекторий', fontsize=14)
plt.legend(loc='best', fontsize=10)
plt.grid(True, linestyle=':', alpha=0.5)
plt.axis('equal')
plt.tight_layout()
plt.show()

# =============================================
# График динамики широты
# =============================================
plt.figure(figsize=(14, 5))
plt.plot(time_ref, lat_ref, 'r-', label='Истинные показания')
plt.plot(time_calc, lat_calc, 'b--', label='Рассчитанные показания')
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Широта [град]', fontsize=12)
plt.title('Сравнение широты', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()

# =============================================
# График динамики долготы
# =============================================
plt.figure(figsize=(14, 5))
plt.plot(time_ref, lon_ref, 'r-', label='Истинные показания')
plt.plot(time_calc, lon_calc, 'b--', label='Рассчитанные показания')
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Долгота [град]', fontsize=12)
plt.title('Сравнение долготы', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(fontsize=10)
plt.tight_layout()
plt.show()

# =============================================
# График VE (столбец 11 vs 7)
# =============================================
plt.figure(figsize=(14, 5))
plt.plot(time_ref, ref_val11, 'r-', label='Истинные показания')
plt.plot(time_calc, calc_val7, 'b--', label='Рассчитанные показания')
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Скорость, м/с', fontsize=12)
plt.title('Сравнение VE', fontsize=14)
plt.legend(loc='upper right', fontsize=10)
plt.grid(True, linestyle=':', alpha=0.5)
plt.ylim(min(ref_val11.min(), calc_val7.min()) - 0.3,
         max(ref_val11.max(), calc_val7.max()) + 0.3)
plt.tight_layout()
plt.show()

# =============================================
# График VN (столбец 12 vs 8)
# =============================================
plt.figure(figsize=(14, 5))
plt.plot(time_ref, ref_val12, 'g-', label='Истинные показания')
plt.plot(time_calc, calc_val8, 'm--', label='Рассчитанные показания')
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Скорость, м/с', fontsize=12)
plt.title('Сравнение VN', fontsize=14)
plt.legend(loc='upper right', fontsize=10)
plt.grid(True, linestyle=':', alpha=0.5)
plt.ylim(min(ref_val12.min(), calc_val8.min()) - 0.3,
         max(ref_val12.max(), calc_val8.max()) + 0.3)
plt.tight_layout()
plt.show()

# =============================================
# График крена с допусками
# =============================================
plt.figure(figsize=(14, 5))
plt.plot(time_ref, roll_ref, 'r-', label='Истинные показания')
plt.plot(time_ref, roll_ref + 0.1, 'k--', linewidth=1, label='Допуск ±0.1°')
plt.plot(time_ref, roll_ref - 0.1, 'k--', linewidth=1)
plt.plot(time_calc, roll_calc, 'b--', label='Рассчитанные показания')
plt.ylabel('Крен [град]', fontsize=12)
plt.title('Сравнение крена', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(fontsize=10)
plt.ylim(roll_ref.min()-0.5, roll_ref.max()+0.5)
plt.tight_layout()
plt.show()

# =============================================
# График тангажа с допусками
# =============================================
plt.figure(figsize=(14, 5))
plt.plot(time_ref, pitch_ref, 'r-', label='Истинные показания')
plt.plot(time_ref, pitch_ref + 0.1, 'k--', linewidth=1, label='Допуск ±0.1°')
plt.plot(time_ref, pitch_ref - 0.1, 'k--', linewidth=1)
plt.plot(time_calc, pitch_calc, 'b--', label='Рассчитанные показания')
plt.ylabel('Тангаж [град]', fontsize=12)
plt.title('Сравнение тангажа', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(fontsize=10)
plt.ylim(pitch_ref.min()-0.5, pitch_ref.max()+0.5)
plt.tight_layout()
plt.show()

# =============================================
# График курса с допусками
# =============================================
plt.figure(figsize=(14, 5))
plt.plot(time_ref, yaw_ref, 'r-', label='Истинные показания')
plt.plot(time_ref, yaw_ref + 1, 'k--', linewidth=1, label='Допуск ±1°')
plt.plot(time_ref, yaw_ref - 1, 'k--', linewidth=1)
plt.plot(time_calc, yaw_calc, 'b--', label='Рассчитанные показания')
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Курс [град]', fontsize=12)
plt.title('Сравнение курса', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(fontsize=10)
plt.ylim(yaw_ref.min()-2, 365)
plt.tight_layout()
plt.show()

# =============================================
# Графики ошибок угловых параметров
# =============================================

# График ошибки крена
plt.figure(figsize=(14, 5))
plt.plot(time_ref, diff_roll, 'r-', label='Ошибка крена')
plt.axhline(y=0.1, color='k', linestyle='--', linewidth=1, label='Допуск ±0.1°')
plt.axhline(y=-0.1, color='k', linestyle='--', linewidth=1)
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Ошибка [град]', fontsize=12)
plt.title('Ошибка крена', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(loc='upper right', fontsize=10)
plt.ylim(-0.15, 0.15)
plt.tight_layout()
plt.show()

# График ошибки тангажа
plt.figure(figsize=(14, 5))
plt.plot(time_ref, diff_pitch, 'b-', label='Ошибка тангажа')
plt.axhline(y=0.1, color='k', linestyle='--', linewidth=1, label='Допуск ±0.1°')
plt.axhline(y=-0.1, color='k', linestyle='--', linewidth=1)
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Ошибка [град]', fontsize=12)
plt.title('Ошибка тангажа', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(loc='upper right', fontsize=10)
plt.ylim(-0.15, 0.15)
plt.tight_layout()
plt.show()

# График ошибки курса
plt.figure(figsize=(14, 5))
plt.plot(time_ref, diff_yaw, 'g-', label='Ошибка курса')
plt.axhline(y=1.0, color='k', linestyle='--', linewidth=1, label='Допуск ±1.0°')
plt.axhline(y=-1.0, color='k', linestyle='--', linewidth=1)
plt.xlabel('Время [с]', fontsize=12)
plt.ylabel('Ошибка [град]', fontsize=12)
plt.title('Ошибка курса', fontsize=14)
plt.grid(True, linestyle=':', alpha=0.5)
plt.legend(loc='upper right', fontsize=10)
plt.ylim(-1.5, 1.5)
plt.tight_layout()
plt.show()

# =============================================
# График высоты (раскомментируйте при необходимости)
# =============================================
# plt.figure(figsize=(14, 5))
# plt.plot(time_ref, alt_ref, 'darkblue', linewidth=1.5)
# plt.xlabel('Время [с]', fontsize=12)
# plt.ylabel('Высота [м]', fontsize=12)
# plt.title('Профиль высоты полёта', fontsize=14)
# plt.grid(True, linestyle=':', alpha=0.5)
# plt.tight_layout()
# plt.show()