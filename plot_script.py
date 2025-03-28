import pandas as pd
import matplotlib.pyplot as plt

def read_column(filename, column_idx):
    """Чтение указанного столбца из файла с обработкой ошибок"""
    try:
        df = pd.read_csv(filename, sep='\s+', header=None, engine='python')
        if df.shape[1] < column_idx + 1:
            print(f"В файле {filename} недостаточно столбцов")
            return None
        return df.iloc[:, column_idx]
    except FileNotFoundError:
        print(f"Файл {filename} не найден")
        return None

# Читаем данные из обоих файлов
col4_file1 = read_column('calculated_data.txt', 3)  # 4-й столбец (индекс 3)
col13_file2 = read_column('20170922_Cessna172_200Hz_Ref.txt', 12)  # 13-й столбец (индекс 12)

# Проверяем успешность чтения данных
if col4_file1 is None or col13_file2 is None:
    exit()

# Создаем подписи для легенды
labels = [
    '4-й столбец из calculated_data.txt',
    '13-й столбец из 20170922_Cessna...'
]

# Создаем график
plt.figure(figsize=(12, 6))

# Рисуем оба графика
plt.plot(col4_file1, marker='', linestyle='-', linewidth=1, color='blue', label=labels[0])
plt.plot(col13_file2, marker='', linestyle='-', linewidth=1, color='red', label=labels[1])

# Настраиваем отображение
plt.title('Сравнение данных из двух файлов', fontsize=14, pad=20)
plt.xlabel('Номер строки', fontsize=12)
plt.ylabel('Значения', fontsize=12)
plt.grid(True, alpha=0.4)
plt.legend(fontsize=10)

# Автоматическая настройка осей для лучшего отображения
plt.tight_layout()
plt.show()