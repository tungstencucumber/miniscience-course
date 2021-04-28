import pandas as pd
from CalcMesh import *

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()

# Размер расчётной сетки, точек на сторону
size = 400
# Шаг точек по пространству
h = 0.001
# Шаг по времени
tau = 0.001
# Количество шагов расчёта
numberOfSteps = 100
# Геометрический центр решётки
center = np.array([0.03, 0.03])
# Рабочая длина волны
la = 0.006
# Скорость распространения волн
speed = 0.5
# Начальная амплитуда
E0 = 1
# Количество излучателей
N = 11
# Расстояние между излучателями
distance = 0.003
# Главное направление решётки
thetaM = np.pi/4

# Создаём сетку заданного размера
MyArray = PhasedArray(center, thetaM, N, distance, E0, la, speed)

MyArray.LinearWave(np.pi/4)
m = CalcMesh(size, h)
# Пишем её начальное состояние в VTK
m.Snapshot(0, "BeamWidthToAzimuth")

printProgressBar(1, numberOfSteps, prefix = 'Progress:', suffix = 'Complete', length = 50)

# Делаем шаги по времени,
# на каждом шаге считаем новое состояние и пишем его в VTK

for i in range(1, numberOfSteps):
    m.DoTimeStep(tau*i, MyArray.GetEmitters())
    m.Snapshot(i, "BeamWidthToAzimuth")
    if i >= int(0.2 * numberOfSteps) and i < int(0.8 * numberOfSteps):
        MyArray.LinearWave(np.pi/4)

    printProgressBar(i + 1, numberOfSteps, prefix = 'Progress:', suffix = 'Complete', length = 50)

emitterArray = np.array([el.Center for el in MyArray.Emitters])
res = np.concatenate((emitterArray, np.array([np.zeros( ( len(MyArray.Emitters) ) )]).T), axis=1 )
df = pd.DataFrame(res)
df.to_csv('LinearWaveEmittersCoordinates.csv')
