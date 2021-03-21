//
// Created by jie zou on 2020/8/7.
//
#include "correspond.h"

int main() {
  const std::string cloth = DATA_PATH"/cloth.obj";
  const std::string body = DATA_PATH"/body.obj";
  const std::string mapfile = DATA_PATH"/cloth2body.map";

  pmp_jie::Correspond correspond;
  correspond.Load(cloth);
  correspond.LoadTarget(body);
  correspond.Deform(mapfile);
  correspond.WriteDeformedMesh();

  return 0;
}
