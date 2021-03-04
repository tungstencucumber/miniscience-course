#include <iostream>
#include <set>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

class CalculableNode // класс расчётной точки
{

friend class CalculableMesh;

protected:
  double x, y, z; // координаты точки
  double vx, vy, vz; // скорость точки
  double density; // плотность жидкости известного рода в данной точке, г/см^3
public:
  // конструктор по умолчанию
  CalculableNode() : x(0.0), y(0.0), z(0.0), vx(0.0), vy(0.0), vz(0.0), density(0.0)
  {}

  // нормальный конструктор со всеми параметрами
  CalculableNode(double x, double y, double z, double vx, double vy, double vz, double smth)
      : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), density(density)
  {}

  void move(double tau) // Метод перемещения точки за время tau с текущей скоростью
  {
    x += vx * tau;
    y += vy * tau;
    z += vz * tau;
  }
};

class Element // класс элемента сетки
{
friend class CalculableMesh;

protected:
  unsigned long NodesIDs[4]; // индексы узлов, составляющих данный элемент сетки
};

class CalculableMesh // класс расчётной сеточки
{
protected:
  std::vector<CalculableNode> MeshNodes; // все узлы сеточки
  std::vector<Element> MeshElements; // все элементы сеточки
public:
  // конструктор сетки по заданной stl-модели
  CalculableMesh(const std::vector<double>& STLNodesCoords, const std::vector<size_t>& STLElementsPoints)
  {
    // загружаем в сеточку точки STL-ной модели
    MeshNodes.resize(STLNodesCoords.size() / 3);
    for(unsigned int i = 0; i < STLNodesCoords.size() / 3; ++i)
    {
      double X = STLNodesCoords[i];
      double Y = STLNodesCoords[i + 1];
      double Z = STLNodesCoords[i + 2];
      double VX = 0.;
      double VY = 0.;
      double VZ = 0.;
      double DENSITY = 0; // тут будет распределение, но попозже
      MeshNodes[i] = CalculableNode(X, Y, Z, VX, VY, VZ, DENSITY);
    }

    // загружаем в сеточку элементы STL-ной модели
    MeshElements.resize(STLElementsPoints.size() / 4);
    for(unsigned int i = 0; i < STLElementsPoints.size() / 4; ++i)
    {
      MeshElements[i].NodesIDs[0] = STLElementsPoints[i*4] - 1;
      MeshElements[i].NodesIDs[1] = STLElementsPoints[i*4 + 1] - 1;
      MeshElements[i].NodesIDs[2] = STLElementsPoints[i*4 + 2] - 1;
      MeshElements[i].NodesIDs[3] = STLElementsPoints[i*4 + 3] - 1;
    }
  }

  // метод изменения сетки во времени с шагом tau
  void DoTimeStep(double tau)
  {
    for(unsigned int i = 0; i < MeshNodes.size(); ++i)
    {
      MeshNodes[i].move(tau);
    }
  }

  // метод записи состояния сетки в формате VTK
  void Snapshot(unsigned int SnapNumber)
  {
    // объявим сетку и её точки
    vtkSmartPointer<vtkUnstructuredGrid> UnstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> DumpPoints = vtkSmartPointer<vtkPoints>::New();

    // зададим скалярное поле плотности жидкости на точках
    auto Density = vtkSmartPointer<vtkDoubleArray>::New();
    Density->SetName("Density");

    // и векторное поле скоростей точек
    auto Velosity = vtkSmartPointer<vtkDoubleArray>::New();
    Velosity->SetName("Velosity");
    Velosity->SetNumberOfComponents(3);

    // выгрузим в VTK-шные точки наши точки
    for(unsigned int i = 0; i < MeshNodes.size(); ++i)
    {
      DumpPoints->InsertNextPoint(MeshNodes[i].x, MeshNodes[i].y, MeshNodes[i].z);
      Density->InsertNextValue(MeshNodes[i].density);
      double _Velosity[3] = {MeshNodes[i].vx, MeshNodes[i].vy, MeshNodes[i].vz};
      Velosity->InsertNextTuple(_Velosity);
    }

    // присоединим VTK-шные точки и оба наших поля в точках к VTK-шной сетке
    UnstructuredGrid->SetPoints(DumpPoints);
    UnstructuredGrid->GetPointData()->AddArray(Velosity);
    UnstructuredGrid->GetPointData()->AddArray(Density);

    // выгружаем в VTK элементы сетки и таким образом объединяем VTK-шные точки в тетраэдры
    for(unsigned int i = 0; i < MeshElements.size(); ++i)
    {
      auto Tetra = vtkSmartPointer<vtkTetra>::New();
      for(unsigned int j = 0; j < 4; ++j)
        Tetra->GetPointIds()->SetId(j, MeshElements[i].NodesIDs[j]);
    }

    // таки сохраним полученную сетку в снапшот
    std::string FileName = "wine-step" + std::to_string(SnapNumber) + ".vtu";
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> Writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    Writer->SetFileName(FileName.c_str());
    Writer->SetInputData(UnstructuredGrid);
    Writer->Write();
  }
};

int main(int argc, char const *argv[]) {
  const unsigned int GMSH_TETR_CODE = 4; // индекс элементов, соответствующий тетраэдрам
  double h = 4.0; // шаг точек по пространству
  double tau = 0.01; // шаг по времени

  gmsh::initialize();
  gmsh::model::add("Wine");

  try {
    gmsh::merge("wine.stl");
  } catch(...) {
    gmsh::logger::write("Could not load STL mesh: bye!");
    gmsh::finalize();
    return -1;
  }

  // Извлекаем геометрию объекта из STL-ины
  double angle = 1;
  bool forceParametrizablePatches = false;
  bool includeBoundary = true;
  double curveAngle = 180;
  gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
  gmsh::model::mesh::createGeometry();

  // Задаем внутренний объём модели
  std::vector<std::pair<int, int> > s;
  gmsh::model::getEntities(s, 2);
  std::vector<int> sl;
  for(auto surf : s) sl.push_back(surf.second);
  int l = gmsh::model::geo::addSurfaceLoop(sl);
  gmsh::model::geo::addVolume({l});

  // Сохранимся
  gmsh::model::geo::synchronize();

  // Задаем мелкость сетки
  int f = gmsh::model::mesh::field::add("MathEval");
  gmsh::model::mesh::field::setString(f, "F", "4");
  gmsh::model::mesh::field::setAsBackgroundMesh(f);

  // Строим сетку
  gmsh::model::mesh::generate(3);

  // Выгружаем данные об узлах сетки: координаты, индекс и ещё какую-то хуйню (каким-то образом параметризованные координаты)
  std::vector<size_t> NodesTags;
  std::vector<double> NodesCoords;
  std::vector<double> ParametricCoords;
  gmsh::model::mesh::getNodes(NodesTags, NodesCoords, ParametricCoords);

  // Выгружаем данные о тетраэдрах
  std::vector<size_t>* TetrasNodesTag = nullptr;
  std::vector<int> ElementsTypes;
  std::vector<std::vector<size_t>> ElementsTags;
  std::vector<std::vector<size_t>> ElementsNodesTags;
  gmsh::model::mesh::getElements(ElementsTypes, ElementsTags, ElementsNodesTags);
  for(unsigned int i = 0; i < ElementsTypes.size(); ++i)
  {
    if(ElementsTypes[i] != GMSH_TETR_CODE)
      continue;
    TetrasNodesTag = &ElementsNodesTags[i];
  }
  // Если мы тетраэдры не нашли, то падаем со свинячим визгом
  if(TetrasNodesTag == nullptr)
  {
    cout << "Could not find tetra data. Aborting calculation." << endl;
    gmsh::finalize();
    return -2;
  }
  cout << "The model has " <<  NodesTags.size() << " nodes and " << TetrasNodesTag->size() / 4 << " tetrs." << endl;

  for(int i = 0; i < NodesTags.size(); ++i)
  {
      assert(i == NodesTags[i] - 1);
  }
  assert(TetrasNodesTag->size() % 4 == 0);
  // Finally...
  CalculableMesh MyMesh(NodesCoords, *TetrasNodesTag);

  gmsh::finalize();

  MyMesh.Snapshot(0);

  return 0;
}
