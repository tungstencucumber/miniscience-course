#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <gmsh.h>

#ifndef PI
#define PI 3.1415927
#endif

int main(int argc, char const *argv[]) {
  gmsh::initialize();
  gmsh::model::add("torus");

  double lc = 1e-1;
  double R = 5;
  double r = 2;
  double r0 = 1.5;
  int count = 0;
  int dim = 10;

  // Outer Ring
  // creating points
  for(int i = 0; i < 2*dim; ++i) {
      for(int j = 0; j < 2*dim; ++j) {
          double psi = -PI + PI*i/dim;
          double phi = PI*j/dim;
          double x = (R + r*cos(psi))*cos(phi);
          double y = (R + r*cos(psi))*sin(phi);
          double z = r*sin(psi);
          gmsh::model::geo::addPoint(x, y, z, lc, j + 2*dim*i + 1);
      }
  }
  // connecting points
  for(int i = 0; i < 2*dim; ++i) {
      for(int j = 0; j < 2*dim - 1; ++j) {
          gmsh::model::geo::addLine(j + 2*dim*i + 1, j + 2*dim*i + 2, j + 2*dim*i + 1);
      }
      gmsh::model::geo::addLine(2*dim*(i + 1), 2*dim*i + 1, 2*dim*(i + 1));
  }
  for(int j = 0; j < 2*dim; ++j) {
      for(int i = 0; i < 2*dim - 1; ++i) {
          gmsh::model::geo::addLine(j + 2*dim*i + 1, j + 2*dim*(i + 1) + 1, i + 2*dim*j + 1 + 4*dim*dim);
      }
      gmsh::model::geo::addLine(j + 2*dim*(2*dim - 1) + 1, j + 1, 2*dim + 2*dim*j + 4*dim*dim);
  }
  // looping edges and making surfaces
  for(int i = 0; i < 2*dim; i++) {
      for(int j = 0; j < 2*dim; j++) {
          gmsh::model::geo::addCurveLoop({(j + 2*dim*i)%(4*dim*dim) + 1, (i + 2*dim*(j + 1))%(4*dim*dim) + 1 + 4*dim*dim, -((j + 2*dim*(i + 1))%(4*dim*dim) + 1), -((i + 2*dim*j)%(4*dim*dim) + 1 + 4*dim*dim)}, j + 2*dim*i + 1);
          gmsh::model::geo::addPlaneSurface({j + 2*dim*i + 1}, j + 2*dim*i + 1);
      }
  }
  std::vector<int> v;
  for(int i = 1; i < 4*dim*dim + 1; i++)
    v.push_back(i);
  gmsh::model::geo::addSurfaceLoop({v}, 1);

  // Inner ring
  // creating points
  for(int i = 0; i < 2*dim; ++i) {
      for(int j = 0; j < 2*dim; ++j) {
          double psi = -PI + PI*i/dim;
          double phi = PI*j/dim;
          double x = (R + r0*cos(psi))*cos(phi);
          double y = (R + r0*cos(psi))*sin(phi);
          double z = r0*sin(psi);
          gmsh::model::geo::addPoint(x, y, z, lc, j + 2*dim*i + 1 + 8*dim*dim);
      }
  }
  // connecting points
  for(int i = 0; i < 2*dim; ++i) {
      for(int j = 0; j < 2*dim - 1; ++j) {
          gmsh::model::geo::addLine(j + 2*dim*i + 1 + 8*dim*dim, j + 2*dim*i + 2 + 8*dim*dim, j + 2*dim*i + 1 + 8*dim*dim);
      }
      gmsh::model::geo::addLine(2*dim*(i + 1) + 8*dim*dim, 2*dim*i + 1 + 8*dim*dim, 2*dim*(i + 1) + 8*dim*dim);
  }
  for(int j = 0; j < 2*dim; ++j) {
      for(int i = 0; i < 2*dim - 1; ++i) {
          gmsh::model::geo::addLine(j + 2*dim*i + 1 + 8*dim*dim, j + 2*dim*(i + 1) + 1 + 8*dim*dim, i + 2*dim*j + 1 + 4*dim*dim + 8*dim*dim);
      }
      gmsh::model::geo::addLine(j + 2*dim*(2*dim - 1) + 1 + 8*dim*dim, j + 1 + 8*dim*dim, 2*dim + 2*dim*j + 4*dim*dim + 8*dim*dim);
  }
  // looping edges and making surfaces
  for(int i = 0; i < 2*dim; i++) {
      for(int j = 0; j < 2*dim; j++) {
          gmsh::model::geo::addCurveLoop({(j + 2*dim*i)%(4*dim*dim) + 1 + 8*dim*dim, (i + 2*dim*(j + 1))%(4*dim*dim) + 1 + 4*dim*dim + 8*dim*dim, -((j + 2*dim*(i + 1))%(4*dim*dim) + 1 + 8*dim*dim), -((i + 2*dim*j)%(4*dim*dim) + 1 + 4*dim*dim + 8*dim*dim)}, j + 2*dim*i + 1 + 8*dim*dim);
          gmsh::model::geo::addPlaneSurface({j + 2*dim*i + 1}, j + 2*dim*i + 1 + 8*dim*dim);
      }
  }
  v.clear();
  for(int i = 1 + 8*dim*dim; i < 4*dim*dim + 1 + 8*dim*dim; i++)
    v.push_back(-i);
  gmsh::model::geo::addSurfaceLoop({v}, 2);
  gmsh::model::geo::addVolume({1, 2}, 1);

  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(3);
  gmsh::write("torus.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}
