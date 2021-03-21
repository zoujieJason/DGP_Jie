//
// Created by jie zou on 2020/8/7.
//
#include "halfedge.h"
pmp_jie::HMesh::HMesh()=default;
int pmp_jie::HMesh::Read(const std::string &file) {
  std::ifstream ifile(file);
  if(ifile.fail()) {
    std::cout << "File read failed. __LINE__ : " << __LINE__ << std::endl;
    return 0;
  }

  std::string buffer;
  while (std::getline(ifile, buffer)) {
    std::istringstream str(buffer);
    std::string symbol;
    str >> symbol;
    if (symbol == "#") {
      continue;
    }
    if (symbol == "v") {
      std::string x, y, z;
      str >> x >> y >> z;
      position_container_.emplace_back(std::stof(x), std::stof(y), std::stof(z));
      auto vertex_ptr = std::make_shared<Vertex>();
      vertex_ptr->idx = NumVerts();
      vertices_container_.push_back(vertex_ptr);
      continue;
    }
    if (symbol == "vn") {
      std::string x, y, z;
      str >> x >> y >> z;
      normal_container_.emplace_back(std::stof(x), std::stof(y), std::stof(z));
      continue;
    }
    if (symbol == "f") {
      std::string v1, v2, v3;
      str >> v1 >> v2 >> v3;
      AddFace(std::stoi(v1) - 1, std::stoi(v2) - 1, std::stoi(v3) - 1);
      continue;
    }
  }
  return 0;
}
int pmp_jie::HMesh::AddFace(const int &v1, const int &v2, const int &v3) {
  const int face_idx = NumFaces();
  auto edge12 = GetEdge(v1, v2);
  auto edge23 = GetEdge(v2, v3);
  auto edge31 = GetEdge(v3, v1);
  assert(!edge12.expired() && !edge23.expired() && !edge31.expired());
  auto halfedge12 = GetHalfedge(edge12, v1, v2);
  auto halfedge23 = GetHalfedge(edge23, v2, v3);
  auto halfedge31 = GetHalfedge(edge31, v3, v1);
  assert(!halfedge12.expired() && !halfedge23.expired() && !halfedge31.expired());
  ConnectTriangle(halfedge12, halfedge23, halfedge31, face_idx);
  vertices_container_[v1]->halfedge = halfedge12;
  vertices_container_[v2]->halfedge = halfedge23;
  vertices_container_[v3]->halfedge = halfedge31;
  AddTriangle(face_idx, v1);
  return 1;
}
size_t pmp_jie::HMesh::NumVerts() const {
  return vertices_container_.size();
}
size_t pmp_jie::HMesh::NumFaces() const {
  return faces_container_.size();
}
bool pmp_jie::HMesh::IsBoundary(const size_t &vid) const {
  assert(vid >=0 && vid < NumVerts());
  const std::weak_ptr<Halfedge> vptr = vertices_container_[vid]->halfedge;
  auto next_ptr = vptr.lock()->opposite_halfedge.lock()->next_halfedge;
  while (next_ptr.lock() != vptr.lock() && next_ptr.lock() != nullptr) {
    next_ptr = next_ptr.lock()->opposite_halfedge.lock()->next_halfedge;
  }
  return !(next_ptr.lock() == vptr.lock());
}
std::vector<std::weak_ptr<pmp_jie::Halfedge>> pmp_jie::HMesh::Halfedges(const int &vid) const {
  std::vector<std::weak_ptr<Halfedge>> result;
  const std::weak_ptr<Halfedge> vptr = vertices_container_[vid]->halfedge;
  result.push_back(vptr);

  auto next_ptr = vptr.lock()->opposite_halfedge.lock()->next_halfedge;
  while (next_ptr.lock() != vptr.lock() && next_ptr.lock() != nullptr) {
    result.push_back(next_ptr);
    next_ptr = next_ptr.lock()->opposite_halfedge.lock()->next_halfedge;
  }

  if (next_ptr.lock() == vptr.lock()) {
    return result;
  }

  auto prev_ptr = vptr.lock()->prev_halfedge;
  while (prev_ptr.lock() != vptr.lock() && prev_ptr.lock() != nullptr) {
    result.push_back(prev_ptr.lock()->opposite_halfedge);
    prev_ptr = prev_ptr.lock()->opposite_halfedge.lock()->prev_halfedge;
  }

  return result;
}
std::vector<std::weak_ptr<pmp_jie::Halfedge>> pmp_jie::HMesh::Verts(const int &fid) const {
  assert(fid >= 0 && fid < NumFaces());
  std::vector<std::weak_ptr<Halfedge>> result;
  const std::weak_ptr<Halfedge> ptr = faces_container_[fid]->halfedge;
  result.push_back(ptr);
  auto next = ptr.lock()->next_halfedge;
  while (next.lock() != ptr.lock()) {
    result.push_back(next);
    next = next.lock()->next_halfedge;
  }
  return result;
}
std::weak_ptr<pmp_jie::Halfedge> pmp_jie::HMesh::GetHalfedge(const std::weak_ptr<Edge> &edge,
                                                             const int &from_vid,
                                                             const int &to_vid) const {
  assert(!edge.expired());
  auto fhalfedge = edge.lock()->from_halfedge;
  auto thalfedge = edge.lock()->to_halfedge;

  assert(!fhalfedge.expired() && !thalfedge.expired());
  if(fhalfedge.lock()->vertex_idx == to_vid) {
    if (thalfedge.lock()->vertex_idx == from_vid) {
      return fhalfedge;
    } else {
      return std::weak_ptr<Halfedge>();
    }
  } else {
    if (thalfedge.lock()->vertex_idx == to_vid && fhalfedge.lock()->vertex_idx == from_vid) {
      return thalfedge;
    } else {
      return std::weak_ptr<Halfedge>();
    }
  }
}
std::weak_ptr<pmp_jie::Edge> pmp_jie::HMesh::FindEdge(const int &from_vid, const int &to_vid) const {
  if (edges_container_map_.count(from_vid)) {
    std::weak_ptr<Edge> ptr = edges_container_map_.find(from_vid)->second;
    while (ptr.lock() != nullptr) {
      if (ptr.lock()->from_halfedge.lock()->vertex_idx == to_vid) {
        return ptr;
      }
      ptr = ptr.lock()->next_edge;
    }
  }
  return std::weak_ptr<Edge>();
}
std::weak_ptr<pmp_jie::Edge> pmp_jie::HMesh::GetEdge(const int &from_vid, const int &to_vid) {
  if (vertices_container_[from_vid]->halfedge.expired() ||
      vertices_container_[to_vid]->halfedge.expired()) {
    return AddEdge(from_vid, to_vid);
  }
  if (!vertices_container_[from_vid]->halfedge.expired() &&
      !vertices_container_[to_vid]->halfedge.expired()) {
    auto from_edge = FindEdge(from_vid, to_vid);
    auto to_edge = FindEdge(to_vid, from_vid);
    assert(from_edge.expired() || to_edge.expired());
    if (from_edge.expired() && !to_edge.expired()) {
      return to_edge;
    } else if (!from_edge.expired() && to_edge.expired()) {
      return from_edge;
    } else {
      return AddEdge(from_vid, to_vid);
    }
  }
  return std::weak_ptr<Edge>();
}
std::weak_ptr<pmp_jie::Edge> pmp_jie::HMesh::AddEdge(const int &v1, const int &v2) {
  auto v1v2_ptr = std::make_shared<Halfedge>();
  auto v2v1_ptr = std::make_shared<Halfedge>();
  v1v2_ptr->vertex_idx = v2;
  v2v1_ptr->vertex_idx = v1;
  v1v2_ptr->opposite_halfedge = v2v1_ptr;
  v2v1_ptr->opposite_halfedge = v1v2_ptr;

  auto edge = std::make_shared<pmp_jie::Edge>();
  edge->from_halfedge = v1v2_ptr;
  edge->to_halfedge = v2v1_ptr;

  v1v2_ptr->edge = edge;
  v2v1_ptr->edge = edge;

  if (!edges_container_map_.count(v1)) {
    edges_container_map_.insert(std::pair<int, std::shared_ptr<Edge>>(v1, edge));
  } else {
    auto ptr = edges_container_map_.find(v1)->second;
    while (ptr->next_edge != nullptr) {
      ptr = ptr->next_edge;
    }
    ptr->next_edge = edge;
  }
  halfedges_container_.push_back(v1v2_ptr);
  halfedges_container_.push_back(v2v1_ptr);

  return edge;
}
std::weak_ptr<pmp_jie::Face> pmp_jie::HMesh::AddTriangle(const int &fid, const int &vid) {
  auto fptr = std::make_shared<Face>();
  fptr->idx = fid;
  fptr->halfedge = vertices_container_[vid]->halfedge;
  faces_container_.push_back(fptr);
  return fptr;
}
int pmp_jie::HMesh::ConnectTriangle(const std::weak_ptr<Halfedge> &halfedge12,
                                    const std::weak_ptr<Halfedge> &halfedge23,
                                    const std::weak_ptr<Halfedge> &halfedge31,
                                    const int &fid) {
  halfedge12.lock()->face_idx = fid;
  halfedge23.lock()->face_idx = fid;
  halfedge31.lock()->face_idx = fid;

  halfedge12.lock()->next_halfedge = halfedge23;
  halfedge23.lock()->prev_halfedge = halfedge12;

  halfedge23.lock()->next_halfedge = halfedge31;
  halfedge31.lock()->prev_halfedge = halfedge23;

  halfedge31.lock()->next_halfedge = halfedge12;
  halfedge12.lock()->prev_halfedge = halfedge31;
  return 1;
}
