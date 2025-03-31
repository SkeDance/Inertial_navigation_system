import numpy as np
import matplotlib.pyplot as plt

# Загрузка расчетных данных
calc_data = np.loadtxt('output.txt', skiprows=2, delimiter='\t')
time_calc = calc_data[:, 0] - calc_data[0, 0]  # Нормирование времени
lat_calc = calc_data[:, 4]                     # Широта
lon_calc = calc_data[:, 5]                     # Долгота
roll_calc = calc_data[:, 2]                    # Крен (столбец 3)
pitch_calc = calc_data[:, 1]                   # Тангаж (столбец 2)
yaw_calc = calc_data[:, 3]                     # Курс (столбец 4)

# Загрузка референсных данных
ref_data = np.loadtxt('20170922_Cessna172_200Hz_Ref.txt', delimiter=None)
time_ref = ref_data[:, 0] - ref_data[0, 0]     # Нормирование времени
lat_ref = ref_data[:, 12]                      # Широта
lon_ref = ref_data[:, 13]                      # Долгота
alt_ref = ref_data[:, 14]                      # Высота
roll_ref = ref_data[:, 7]                      # Крен
pitch_ref = ref_data[:, 8]                     # Тангаж
yaw_ref = ref_data[:, 9]                       # Курс

# =============================================
# Графики траекторий
# =============================================
plt.figure(figsize=(12, 10))

# Траектория в координатах
plt.subplot(3, 1, 1)
plt.plot(lon_ref, lat_ref, 'r-', label='Референс')
plt.plot(lon_calc, lat_calc, 'b--', label='Расчёт')
plt.xlabel('Долгота [град]')
plt.ylabel('Широта [град]')
plt.title('Сравнение траекторий')
plt.legend()
plt.grid(True)

# Широта во времени
plt.subplot(3, 1, 2)
plt.plot(time_ref, lat_ref, 'r-', label='Референс')
plt.plot(time_calc, lat_calc, 'b--', label='Расчёт')
plt.xlabel('Время [с]')
plt.ylabel('Широта [град]')
plt.title('Изменение широты')
plt.legend()
plt.grid(True)

# Долгота во времени
plt.subplot(3, 1, 3)
plt.plot(time_ref, lon_ref, 'r-', label='Референс')
plt.plot(time_calc, lon_calc, 'b--', label='Расчёт')
plt.xlabel('Время [с]')
plt.ylabel('Долгота [град]')
plt.title('Изменение долготы')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# =============================================
# Графики ориентации (сравнение)
# =============================================
plt.figure(figsize=(12, 10))

# Крен
plt.subplot(3, 1, 1)
plt.plot(time_ref, roll_ref, 'r-', label='Референс')
plt.plot(time_calc, roll_calc, 'b--', label='Расчёт')
plt.ylabel('Крен [град]')
plt.title('Сравнение углов ориентации')
plt.legend()
plt.grid(True)

# Тангаж
plt.subplot(3, 1, 2)
plt.plot(time_ref, pitch_ref, 'r-', label='Референс')
plt.plot(time_calc, pitch_calc, 'b--', label='Расчёт')
plt.ylabel('Тангаж [град]')
plt.legend()
plt.grid(True)

# Курс
plt.subplot(3, 1, 3)
plt.plot(time_ref, yaw_ref, 'r-', label='Референс')
plt.plot(time_calc, yaw_calc, 'b--', label='Расчёт')
plt.xlabel('Время [с]')
plt.ylabel('Курс [град]')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# =============================================
# График высоты
# =============================================
plt.figure(figsize=(10, 5))
plt.plot(time_ref, alt_ref, 'k-')
plt.xlabel('Время [с]')
plt.ylabel('Высота [м]')
plt.title('Профиль высоты')
plt.grid(True)
plt.tight_layout()
plt.show()