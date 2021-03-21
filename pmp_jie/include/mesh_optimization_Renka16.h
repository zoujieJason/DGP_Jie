//
// Created by jie zou on 2021/3/21.
//

#ifndef DGP_JIE_PMP_JIE_INCLUDE_MESH_OPTIMIZATION_RENKA16_H_
#define DGP_JIE_PMP_JIE_INCLUDE_MESH_OPTIMIZATION_RENKA16_H_

#include "polygon_mesh.h"

namespace pmp_jie {

double optimize_vert_angle_based_Renka16(pmp_jie::Base *mesh_ptr, const size_t &vid);

double area_energy(const cinolib::Polygonmesh<> &trimesh, const size_t &vid);

double optimize_vert_uniform_area_Renka16(pmp_jie::Base *mesh_ptr, const size_t &vid);


}
#endif //DGP_JIE_PMP_JIE_INCLUDE_MESH_OPTIMIZATION_RENKA16_H_
