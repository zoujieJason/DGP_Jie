//
// Created by jie zou on 2019/11/7.
//

#include "../include/QEM_structures.h"

pmp_jie::Vertex::Vertex(const Eigen::Vector3f &pos, const int &idRef) {
    position = pos;
    id = idRef;
    cost = 1000000.0f;
}

int pmp_jie::Vertex::removeTriangle(std::weak_ptr<Triangle> tri) {
    assert(!tri.expired());
    for(auto it = face.begin(); it != face.end(); ++it) {
        if((*it).lock() == tri.lock()) {
            face.erase(it);
            break;
        }
    }
    return 1;
}

int pmp_jie::Vertex::removeNeighbor(std::weak_ptr<Vertex> vertex) {
    assert(!vertex.expired());
    for(auto it = neighbor.begin(); it != neighbor.end(); ++it) {
        if((*it).lock() == vertex.lock()) {
            neighbor.erase(it);
            break;
        }
    }
    return 1;
}

int pmp_jie::Vertex::addTriangle(std::weak_ptr<Triangle> tri) {
    assert(!tri.expired());
    if(!hasTriangle(tri)) {
        face.push_back(tri);
    }
    return 0;
}

int pmp_jie::Vertex::addNeighbor(std::weak_ptr<Vertex> vertex) {
    assert(!vertex.expired());
    if(!hasNeighbor(vertex)) {
        neighbor.push_back(vertex);
    }
    return 1;
}

bool pmp_jie::Vertex::hasNeighbor(std::weak_ptr<Vertex> vertex) const {
    assert(!vertex.expired());
    for(auto it = neighbor.begin(); it != neighbor.end(); ++it) {
        if((*it).lock() == vertex.lock()) {
            return true;
        }
    }
    return false;
}


bool pmp_jie::Vertex::hasTriangle(std::weak_ptr<Triangle> tri) const {
    assert(!tri.expired());
    for(auto it = face.begin(); it != face.end(); ++it) {
        if((*it).lock() == tri.lock()) {
            return true;
        }
    }
    return false;
}

Eigen::Matrix4f pmp_jie::Vertex::getMatQ() const {
    Eigen::Matrix4f Q;
    Q.setZero();
    const auto va = position;
    for(const auto &tri: face) {
        const Eigen::Vector3f vn = tri.lock()->normal;
        const float D = vn(0)*va(0) + vn(1)*va(1) + vn(2)*va(2);
        Eigen::Vector4f p(vn(0), vn(1), vn(2), -D);
        p /= vn.norm();
        Q = Q + p * p.transpose();
    }
    return Q;
}

std::vector<std::weak_ptr<pmp_jie::Triangle>>
pmp_jie::Vertex::getTrianglesBetweenVlAndVr(std::weak_ptr<Vertex> vl, std::weak_ptr<Vertex> vr) const {
    assert(!vl.expired() || !vr.expired());
    std::vector<std::weak_ptr<Triangle>> tris;
    bool flag = true;
    while(!vl.expired()) {
        auto tri = getCwTriangleFromVertex(vl);
        if(tri.expired()) break;
        vl = tri.lock()->getPervVertex(vl);
        tris.push_back(tri);
        if(!vr.expired()) {
            if(vl.lock() == vr.lock()) {
                flag = false;
                break;
            }
        }
    }
    if(flag) {
        while(!vr.expired()) {
            auto tri = getCCwTriangleFromVertex(vr);
            if(tri.expired()) break;
            vr = tri.lock()->getNextVertex(vr);
            tris.push_back(tri);
        }
    }
    return tris;
}

std::weak_ptr<pmp_jie::Triangle> pmp_jie::Vertex::getCwTriangleFromVertex(std::weak_ptr<Vertex> vertex) const {
    assert(!vertex.expired());
    std::weak_ptr<Triangle> tri;
    for(const auto &tri: face) {
        if(!tri.lock()->hasVertex(vertex)) continue;
        if(tri.lock()->getNextVertex(vertex).lock()->id == id) {
            return tri;
        }
    }
    return tri;
}

std::weak_ptr<pmp_jie::Triangle> pmp_jie::Vertex::getCCwTriangleFromVertex(std::weak_ptr<Vertex> vertex) const {
    assert(!vertex.expired());
    std::weak_ptr<Triangle> tri;
    for(const auto &tri: face) {
        if(!tri.lock()->hasVertex(vertex)) continue;
        if(tri.lock()->getPervVertex(vertex).lock()->id == id) {
            return tri;
        }
    }
    return tri;
}

int pmp_jie::Vertex::deleteNeighbor(std::weak_ptr<Vertex> vertex) {
    assert(!vertex.expired());
    for(auto it = neighbor.begin(); it != neighbor.end(); ++it) {
        if((*it).lock() == vertex.lock()) {
            neighbor.erase(it);
            break;
        }
    }
    return 0;
}

int pmp_jie::Vertex::deleteTriangle(std::weak_ptr<Triangle> tri) {
    assert(!tri.expired());
    for(auto it = face.begin(); it != face.end(); ++it) {
        if((*it).lock() == tri.lock()) {
            face.erase(it);
            break;
        }
    }
    return 0;
}


pmp_jie::Triangle::Triangle(std::weak_ptr<Vertex> v0Ref,
                   std::weak_ptr<Vertex> v1Ref,
                   std::weak_ptr<Vertex> v2Ref) {
    v0 = v0Ref;
    v1 = v1Ref;
    v2 = v2Ref;
    computeNormal();
}

void pmp_jie::Triangle::computeNormal() {
    const Eigen::Vector3f va(v0.lock()->position);
    const Eigen::Vector3f vb(v1.lock()->position);
    const Eigen::Vector3f vc(v2.lock()->position);
    normal = (vb - va).cross(vc - va);
}

bool pmp_jie::Triangle::hasVertex(std::weak_ptr<Vertex> vertex) const {
    assert(!vertex.expired());
    return (vertex.lock() == v0.lock() || vertex.lock() == v1.lock() || vertex.lock() == v2.lock());
}

int pmp_jie::Triangle::replaceVertex(std::weak_ptr<Vertex> vold, std::weak_ptr<Vertex> vnew) {
    assert(!vold.expired() && !vnew.expired());
    if(vold.lock() == v0.lock()) {
        v0 = vnew;
    }
    else if (vold.lock() == v1.lock()) {
        v1 = vnew;
    }
    else {
        assert(vold.lock() == v2.lock());
        v2 = vnew;
    }
    computeNormal();
    return 0;
}

std::weak_ptr<pmp_jie::Vertex> pmp_jie::Triangle::getNextVertex(std::weak_ptr<Vertex> vertex) const {
    assert(!vertex.expired());
    if(vertex.lock() == v0.lock()) {
        return v1;
    }
    else if(vertex.lock() == v1.lock()) {
        return v2;
    }
    else {
        return v0;
    }
}

std::weak_ptr<pmp_jie::Vertex> pmp_jie::Triangle::getPervVertex(std::weak_ptr<Vertex> vertex) const {
    assert(!vertex.expired());
    if(vertex.lock() == v0.lock()) {
        return v2;
    }
    else if(vertex.lock() == v1.lock()) {
        return v0;
    }
    else {
        return v1;
    }
}

pmp_jie::VSplit::VSplit() {

}

void pmp_jie::VSplit::setVfPos(const Eigen::Vector3f &pos) {
    vfPos = pos;
}

void pmp_jie::VSplit::setVtPos(const Eigen::Vector3f &pos) {
    vtPos = pos;
}

void pmp_jie::VSplit::setVfId(const int &id) {
    indices(0) = id;
}

void pmp_jie::VSplit::setVtId(const int &id) {
    indices(1) = id;
}

void pmp_jie::VSplit::setVlId(const int &id) {
    indices(2) = id;
}

void pmp_jie::VSplit::setVrId(const int &id) {
  indices(3) = id;
}