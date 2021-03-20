//
// Created by jie zou on 2019/11/7.
//

#include "../include/mesh_QEM_simplification.h"

int pmp_jie::Mesh::read(const char *filePath) {
    std::ifstream file(filePath);
    if (file.fail())
        return __LINE__;

    // TODO: What if this file contains ilegal lines?
    std::string buffer;
    int vId = 0;
    while(getline(file, buffer)) {
        std::istringstream strIn(buffer);
        std::string symbol;
        strIn >> symbol;
        if(symbol == "#") {
            continue;
        }
        if(symbol == "v") {
            std::string strData1, strData2, strData3;
            strIn >> strData1 >> strData2 >> strData3;
            auto ptr = std::make_shared<Vertex>(Eigen::Vector3f(std::stof(strData1),
                                                                std::stof(strData2),
                                                                std::stof(strData3)), vId);
            verticesMap.insert(std::pair<int, std::shared_ptr<Vertex>>(vId, ptr));
            permutation.push_back(vId);
            ++vId;
            continue;
        }
        if(symbol == "f") {
            std::string strData1, strData2, strData3;
            strIn >> strData1 >> strData2 >> strData3;
            std::weak_ptr<Vertex> weakPtr1 = verticesMap.find(std::stoi(strData1) - 1)->second;
            std::weak_ptr<Vertex> weakPtr2 = verticesMap.find(std::stoi(strData2) - 1)->second;
            std::weak_ptr<Vertex> weakPtr3 = verticesMap.find(std::stoi(strData3) - 1)->second;
            auto triPtr = std::make_shared<Triangle>(weakPtr1, weakPtr2, weakPtr3);
            triangleSet.insert(triPtr);
            std::weak_ptr<Triangle> triWeakPtr = triPtr;
            verticesMap.find(std::stoi(strData1) - 1)->second->addNeighbor(weakPtr2);
            verticesMap.find(std::stoi(strData1) - 1)->second->addNeighbor(weakPtr3);
            verticesMap.find(std::stoi(strData1) - 1)->second->addTriangle(triWeakPtr);
            verticesMap.find(std::stoi(strData2) - 1)->second->addNeighbor(weakPtr1);
            verticesMap.find(std::stoi(strData2) - 1)->second->addNeighbor(weakPtr3);
            verticesMap.find(std::stoi(strData2) - 1)->second->addTriangle(triWeakPtr);
            verticesMap.find(std::stoi(strData3) - 1)->second->addNeighbor(weakPtr1);
            verticesMap.find(std::stoi(strData3) - 1)->second->addNeighbor(weakPtr2);
            verticesMap.find(std::stoi(strData3) - 1)->second->addTriangle(triWeakPtr);
            continue;
        }
        if(symbol == "vt") {
            continue;
        }
        if(symbol == "vn") {
            continue;
        }
    }
    file.close();
    computeAllCollapse();
    return 0;
}

int pmp_jie::Mesh::simplify(const int &resultNumber) {
    std::cout << "origin mesh vertice number: " << verticesMap.size() << std::endl;
    std::cout << "targer mesh vertice number: " << resultNumber << std::endl;
    std::cout << std::endl;
    auto currentNum = verticesMap.size();
    while(resultNumber != currentNum) {
        edgeCollapse();
        currentNum = verticesMap.size();
        std::cout << currentNum << std::endl;
    }
    return 0;
}

int pmp_jie::Mesh::write(const char *filePath) const {

    std::ofstream file(filePath);
    file.clear();

//    //vertices
//    for(size_t vId = 0; vId < verticesMap.size(); ++vId) {
//        int curCorrespondence = permutation[vId];
//        while(!verticesMap.count(curCorrespondence)) {
//            curCorrespondence = permutation[curCorrespondence];
//        }
//        auto vertex = verticesMap.find(curCorrespondence)->second;
//        vertex->id = vId;
//        file << "v"
//             << " " << vertex->position(0)
//             << " " << vertex->position(1)
//             << " " << vertex->position(2) << std::endl;
//    }

    //vertices
    int vId = 0;
    for(auto &pair: verticesMap) {
        auto vertex = pair.second;
        vertex->id = vId;
        ++vId;
        file << "v"
             << " " << vertex->position(0)
             << " " << vertex->position(1)
             << " " << vertex->position(2) << std::endl;
    }
    //triangles
    for(const auto &triangle: triangleSet) {
        file << "f"
             << " " << (*triangle).v0.lock()->id + 1
             << " " << (*triangle).v1.lock()->id + 1
             << " " << (*triangle).v2.lock()->id + 1 << std::endl;
    }
    //reset
    for(const auto &pair: verticesMap) {
        pair.second->id = pair.first;
    }
    if (file.fail())
        return __LINE__;
    return 1;
}

int pmp_jie::Mesh::computeAllCollapse() {
    for(const auto &pair: verticesMap) {
        computeCollapseAtVertex(pair.second);
    }
    return 0;
}

std::weak_ptr<pmp_jie::Vertex> pmp_jie::Mesh::getMinimumCostVertex() const {
    std::weak_ptr<Vertex> vertex(verticesMap.begin()->second);
    for(const auto &pair: verticesMap) {
        if(pair.second->cost < vertex.lock()->cost) {
            vertex = pair.second;
        }
    }
    return vertex;
}

int pmp_jie::Mesh::edgeCollapse() {
    auto minVertex = getMinimumCostVertex();
    assert(!minVertex.expired());
    assert(!minVertex.lock()->collapseVertex.expired());
    auto minCollapseVertex = minVertex.lock()->collapseVertex;
    setVsplit(minVertex);
    std::cout << "Edge : " << minVertex.lock()->id << " " << minCollapseVertex.lock()->id << std::endl;
    for(const auto &tri: minVertex.lock()->face) {
        std::cout << tri.lock()->v0.lock()->id << " " << tri.lock()->v1.lock()->id << " " << tri.lock()->v2.lock()->id << " || ";
    }
    std::cout << std::endl;
    for(const auto &tri: minCollapseVertex.lock()->face) {
        std::cout << tri.lock()->v0.lock()->id << " " << tri.lock()->v1.lock()->id << " " << tri.lock()->v2.lock()->id << " || ";
    }
    std::cout << std::endl;
    minCollapseVertex.lock()->position = (minVertex.lock()->position + minCollapseVertex.lock()->position) / 2;
    std::vector<std::weak_ptr<Triangle>> trisHasCollapseVertex;
    std::vector<std::weak_ptr<Triangle>> trisNoCollapseVertex;
    for(const auto &tri: minVertex.lock()->face) {
        if (tri.lock()->hasVertex(minCollapseVertex)) {
            trisHasCollapseVertex.push_back(tri);
        }
        else {
            trisNoCollapseVertex.push_back(tri);
        }
    }
    for(const auto &tri: trisNoCollapseVertex) {
        tri.lock()->replaceVertex(minVertex, minCollapseVertex);
        minCollapseVertex.lock()->addTriangle(tri);
        auto v1 = tri.lock()->getNextVertex(minCollapseVertex);
        auto v2 = tri.lock()->getPervVertex(minCollapseVertex);
        v1.lock()->deleteNeighbor(minVertex);
        v2.lock()->deleteNeighbor(minVertex);
        minCollapseVertex.lock()->addNeighbor(v1);
        v1.lock()->addNeighbor(minCollapseVertex);
        minCollapseVertex.lock()->addNeighbor(v2);
        v2.lock()->addNeighbor(minCollapseVertex);
    }
    permutation[minVertex.lock()->id] = verticesMap.size() - 1;
    minCollapseVertex.lock()->deleteNeighbor(minVertex);
    deleteVertex(minVertex);
    for(const auto &tri: trisHasCollapseVertex) {
        deleteTriangle(tri);
    }
    updateCollapseAroundVertex(minCollapseVertex);
    return 0;
}

int pmp_jie::Mesh::computeCollapseAtVertex(std::weak_ptr<Vertex> vertex) {
    assert(!vertex.expired());
    vertex.lock()->cost = 10000000.0f;
    vertex.lock()->collapseVertex = *vertex.lock()->neighbor.begin();
    const auto Q1 = vertex.lock()->getMatQ();
    const auto va = vertex.lock()->position;
    for(const auto &adjacentVertex: vertex.lock()->neighbor) {
        const auto Q = Q1 + adjacentVertex.lock()->getMatQ();
        const auto vmid3f = (va + adjacentVertex.lock()->position) / 2;
        const auto vmid4f = Eigen::Vector4f(vmid3f(0), vmid3f(1), vmid3f(2), 1.0);
        float minCost = vmid4f.transpose() * Q * vmid4f;
        if(vertex.lock()->cost > minCost) {
            vertex.lock()->cost = minCost;
            vertex.lock()->collapseVertex = adjacentVertex;
        }
    }
    return 0;
}

int pmp_jie::Mesh::deleteVertex(std::weak_ptr<Vertex> vertex) {
    assert(!vertex.expired());
    for(const auto &nei: vertex.lock()->neighbor) {
        nei.lock()->removeNeighbor(vertex);
    }
    assert(verticesMap.count(vertex.lock()->id));
    verticesMap.erase(verticesMap.find(vertex.lock()->id));
    return 0;
}

int pmp_jie::Mesh::deleteTriangle(std::weak_ptr<Triangle> tri) {
    assert(!tri.expired());
    if(!tri.lock()->v0.expired()) {
        tri.lock()->v0.lock()->removeTriangle(tri);
    }
    if(!tri.lock()->v1.expired()) {
        tri.lock()->v1.lock()->removeTriangle(tri);
    }
    if(!tri.lock()->v2.expired()) {
        tri.lock()->v2.lock()->removeTriangle(tri);
    }
    if(triangleSet.count(tri.lock())) {
        triangleSet.erase(triangleSet.find(tri.lock()));
    }
    return 0;
}

int pmp_jie::Mesh::updateCollapseAroundVertex(std::weak_ptr<Vertex> vertex) {
    computeCollapseAtVertex(vertex);
    for(const auto &adjacentVertex: vertex.lock()->neighbor) {
        computeCollapseAtVertex(adjacentVertex);
    }
    return 0;
}

int pmp_jie::Mesh::vertexSplit() {
    VSplit vsplit = *vsplits.rbegin();
    vsplits.pop_back();
    std::cout << vsplit.indices.transpose() << std::endl;
    permutation[vsplit.indices(0)] = vsplit.indices(0);
    std::weak_ptr<Vertex> weakVtPtr;
    std::weak_ptr<Vertex> weakVlPtr;
    std::weak_ptr<Vertex> weakVrPtr;

    assert(verticesMap.count(vsplit.indices(1)));
    weakVtPtr = verticesMap.find(vsplit.indices(1))->second;
    weakVtPtr.lock()->position = vsplit.vtPos;
    if(verticesMap.count(vsplit.indices(2))) {
        weakVlPtr = verticesMap.find(vsplit.indices(2))->second;
    }
    if(verticesMap.count(vsplit.indices(3))) {
        weakVrPtr = verticesMap.find(vsplit.indices(3))->second;
    }
    assert(!weakVlPtr.expired() || !weakVrPtr.expired());

    auto newVertexPtr = std::make_shared<Vertex>(vsplit.vfPos, vsplit.indices(0));
    verticesMap.insert(std::pair<int, std::shared_ptr<Vertex>>(vsplit.indices(0), newVertexPtr));
    std::weak_ptr<Vertex> weakNewVertexPtr = newVertexPtr;
    std::vector<std::weak_ptr<Triangle>> tris = weakVtPtr.lock()->getTrianglesBetweenVlAndVr(weakVlPtr, weakVrPtr);
    for(const auto &tri: tris) {
        assert(!tri.expired());
        tri.lock()->replaceVertex(weakVtPtr, weakNewVertexPtr);
        weakNewVertexPtr.lock()->addTriangle(tri);
        weakVtPtr.lock()->deleteTriangle(tri);

        auto v1 = tri.lock()->getNextVertex(weakNewVertexPtr);
        auto v2 = tri.lock()->getPervVertex(weakNewVertexPtr);
        weakNewVertexPtr.lock()->addNeighbor(v1);
        v1.lock()->addNeighbor(weakNewVertexPtr);
        weakNewVertexPtr.lock()->addNeighbor(v2);
        v2.lock()->addNeighbor(weakNewVertexPtr);
    }
    if(!weakVlPtr.expired()) {
        addTriangle(weakNewVertexPtr, weakVlPtr, weakVtPtr);
    }
    if(!weakVrPtr.expired()) {
        addTriangle(weakVtPtr, weakVrPtr, weakNewVertexPtr);
    }
    return 0;
}

int pmp_jie::Mesh::addTriangle(std::weak_ptr<Vertex> v1, std::weak_ptr<Vertex> v2, std::weak_ptr<Vertex> v3) {
    auto triPtr = std::make_shared<Triangle>(v1, v2, v3);
    std::weak_ptr<Triangle> weakTriPtr = triPtr;
    v1.lock()->addTriangle(weakTriPtr);
    v2.lock()->addTriangle(weakTriPtr);
    v3.lock()->addTriangle(weakTriPtr);

    v1.lock()->addNeighbor(v2);
    v1.lock()->addNeighbor(v3);
    v2.lock()->addNeighbor(v1);
    v2.lock()->addNeighbor(v3);
    v3.lock()->addNeighbor(v1);
    v3.lock()->addNeighbor(v2);
    triangleSet.insert(triPtr);
    return 0;
}

int pmp_jie::Mesh::setVsplit(std::weak_ptr<Vertex> vertex) {
    assert(!vertex.expired());
    assert(!vertex.lock()->collapseVertex.expired());
    auto collapseVertex = vertex.lock()->collapseVertex;
    VSplit vsplit;
    vsplit.setVfPos(vertex.lock()->position);
    vsplit.setVtPos(collapseVertex.lock()->position);
    vsplit.setVfId(vertex.lock()->id);
    vsplit.setVtId(collapseVertex.lock()->id);
    int vlId = -1, vrId = -1;
    for(const auto &vertexTri: vertex.lock()->face) {
        for(const auto &collapseVertexTri: collapseVertex.lock()->face) {
            if(vertexTri.lock() == collapseVertexTri.lock()) {
                auto nextVertex = vertexTri.lock()->getNextVertex(vertex);
                if(nextVertex.lock() != collapseVertex.lock()) {
                    vlId = nextVertex.lock()->id;
                }
                else {
                    vrId = (vertexTri.lock()->getNextVertex(collapseVertex)).lock()->id;
                }
            }
        }
    }
    vsplit.setVlId(vlId);
    vsplit.setVrId(vrId);
    vsplits.push_back(vsplit);
    return 0;
}

int pmp_jie::Mesh::refine(const int &resultNumber) {
    std::cout << "current mesh vertice number: " << verticesMap.size() << std::endl;
    std::cout << "targer mesh vertice number: " << resultNumber << std::endl;
    std::cout << std::endl;
    auto currentNum = verticesMap.size();
    while(resultNumber != currentNum) {
        vertexSplit();
        currentNum = verticesMap.size();
        std::cout << "current mesh vertice number: " << currentNum << std::endl;
        std::cout << "remain vsplit number: " << vsplits.size() << std::endl;
    }
    return 0;
}

void pmp_jie::Mesh::printPermutation() const {
    std::cout << "permutation: ";
    for(const auto &val: permutation) {
        std::cout <<  val << " ";
    }
    std::cout << std::endl;
}

int pmp_jie::Mesh::nVertices() const {
    return verticesMap.size();
}
