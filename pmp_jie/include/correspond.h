//
// Created by jie zou on 2020/8/7.
//

#ifndef DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_CORRESPOND_H_
#define DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_CORRESPOND_H_

#include "mesh_io.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include "linear_solver.h"
#include "tools.h"

//CGAL
typedef CGAL::Simple_cartesian<double> K;
typedef typename K::FT FT;
typedef typename K::Ray_3 Ray;
typedef typename K::Point_3 Point_3;
typedef typename K::Triangle_3 Triangle;
typedef typename std::vector<Triangle>::iterator Iterator;
typedef typename std::vector<Triangle>::const_iterator CIterator;
typedef typename CGAL::AABB_triangle_primitive<K, Iterator> Primitive;
typedef typename CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef typename CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

namespace pmp_jie {
class Correspond: public pmp_jie::CTriMesh {
 public:
  Correspond();

  int LoadTarget(const std::string &target);
  int Deform(const std::string &mapfile);
  int WriteCorrespondMesh(const std::string &filename_with_path= "");

 protected:
  int UsrAllocateVars() override;
  int BuildTree();
  int LoadFixedPoints(const std::string &mapfile);

  size_t SetSmoothnessTerm(std::vector<Eigen::Triplet<double>> &triplets) const;
  size_t SetIdentityTerm(std::vector<Eigen::Triplet<double>> &triplets) const;
  int SetIdentityMatrix(const size_t &rows, Eigen::MatrixXd &matrix) const;
  size_t SetConstraintTerm(std::vector<Eigen::Triplet<double>> &triplets) const;
  int SetConstraintMatrix(Eigen::MatrixXd &matrix) const;
  size_t SetPolyTriplets(const size_t &pid,const size_t &curr_row,const double &weight,std::vector<Eigen::Triplet<double>> &triplets) const;
  size_t SetVaildPointsTerm(std::vector<Eigen::Triplet<double>> &triplets) const;

  int UpdateClosestVaildPoints(const Eigen::MatrixXd &curr_result);
  int FindCorrespondence();

  Point_3 GetClosestPoint(const Point_3 &query) const;
  size_t GetClosestPolyId(const Point_3 &query) const;
  int GetClosestVertId(const size_t &cor_pid, const Eigen::Vector3d &pos) const;
  Point_3 Eigen2CGAL(const Eigen::Vector3d &vec) const;
  Eigen::Vector3d CGAL2Eigen(const Point_3 &vec) const;

  cinolib::Trimesh<> target_;
  std::vector<Triangle> tris_;
  Tree tree_;
  std::map<int,int> fixed_points_;
  Eigen::MatrixXd vaild_points_;

  //map
  std::map<size_t, size_t> v2v_;
  std::map<size_t, size_t> p2p_;
  std::map<size_t, double> e2d_;

};
}

#endif //DIGITALGEOMETRYPROCESSING_DGP_INCLUDE_CORRESPOND_H_
