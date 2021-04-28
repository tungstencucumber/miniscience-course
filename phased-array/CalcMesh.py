import numpy as np
import vtk

class WaveEmitter:

    def __init__(self, c, v, delta0, E0, la):
        self.Center = c
        self.Speed = v
        self.InitialDelta = delta0
        self.InitialAmplitude = E0
        self.Wavelength = la

    def SetDelta(self, delta):
        self.InitialDelta = delta

class PhasedArray:

    def __init__(self, c, thetaM, N, d, E0, la, v):
        self.EmitterDistance = d
        self.EmitterNumber = N
        self.MainDirection = thetaM
        self.ArrayCenter = c
        self.Emitters = []
        for i in range(N):
            self.Emitters.append(WaveEmitter(np.array([c[0] - ((1 - N)/2 + i)*d*np.sin(thetaM), c[1] + ((1 - N)/2 + i)*d*np.cos(thetaM)]), v, 0., E0, la))

    def LinearWave(self, theta):
        for i in range(self.EmitterNumber):
            self.Emitters[i].SetDelta(-i * self.EmitterDistance * np.sin(theta - self.MainDirection))

    def FocusedWave(self, focus):
        ro = np.sqrt((focus[0] - self.ArrayCenter[0])**2 + (focus[1] - self.ArrayCenter[1])**2)
        for i in range(self.EmitterNumber):
            xi = self.ArrayCenter[0] - ((1 - self.EmitterNumber)/2 + i)*self.EmitterDistance*np.sin(self.MainDirection)
            yi = self.ArrayCenter[1] + ((1 - self.EmitterNumber)/2 + i)*self.EmitterDistance*np.cos(self.MainDirection)
            ri = np.sqrt((focus[0] - xi)**2 + (focus[1] - yi)**2)
            self.Emitters[i].SetDelta(ri - ro)

    def GetEmitters(self):
        return self.Emitters

class CalcMesh:

    def __init__(self, size, step):
        self.Nodes = np.mgrid[0:size-1:np.complex(size), 0:size-1:np.complex(size)]
        self.Nodes *= step
        self.Nodes = np.append(self.Nodes, [np.zeros(shape=(size, size), dtype=np.double)], 0)
        self.Size = size
        self.ElectricField = np.zeros(shape=(size, size), dtype=np.double)

    def CalcWave(self, r, t, emitter):
        E = np.zeros(np.shape(self.ElectricField), dtype = np.double)
        x = r[0] - emitter.Center[0]
        y = r[1] - emitter.Center[1]
        E = emitter.InitialAmplitude * np.cos( 2*np.pi/emitter.Wavelength * (emitter.Speed * t - np.sqrt(x**2 + y**2) + emitter.InitialDelta))
        return E

    def DoTimeStep(self, t, emitters):
        self.ElectricField = np.zeros(shape=(self.Size, self.Size), dtype=np.double)
        for Emitter in emitters:
            self.ElectricField += self.CalcWave(self.Nodes, t, Emitter)


    def Snapshot(self, snap_number, name):
        # Сетка в терминах VTK
        structuredGrid = vtk.vtkStructuredGrid()
        # Точки сетки в терминах VTK
        points = vtk.vtkPoints()

        # Векторное поле на точках сетки
        ElectricField = vtk.vtkDoubleArray()
        ElectricField.SetName("ElectricField")

        # Обходим все точки нашей расчётной сетки
        # Делаем это максимально неэффективным, зато наглядным образом
        number = len(self.Nodes[0])
        for i in range(0, number):
            for j in range(0, number):
                points.InsertNextPoint(self.Nodes[0][i,j], self.Nodes[1][i,j], self.Nodes[2][i,j])
                ElectricField.InsertNextValue(self.ElectricField[i,j])

        # Задаём размеры VTK-сетки (в точках, по трём осям)
        structuredGrid.SetDimensions(number, number, 1)
        # Грузим точки в сетку
        structuredGrid.SetPoints(points)

        structuredGrid.GetPointData().AddArray(ElectricField)

        # Создаём снапшот в файле с заданным именем
        writer = vtk.vtkXMLStructuredGridWriter()
        writer.SetInputDataObject(structuredGrid)
        writer.SetFileName(name + "/" + name + "-step-" + str(snap_number) + ".vts")
        writer.Write()
