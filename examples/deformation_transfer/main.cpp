//
// Created by jie zou on 2020/8/7.
//
#include "deformation_transfer.h"

int main() {
  const std::string cloth      = DATA_PATH"/cloth.obj";
  const std::string body       = DATA_PATH"/dress_body.obj";
  const std::string mapfile    = DATA_PATH"/cloth2body.map";

  const std::string pose_body = DATA_PATH"/pose_body.obj";

  pmp_jie::DeformationTransfer transfer;

  //correspond
  transfer.Load(cloth);
  transfer.LoadTarget(body);
  transfer.Deform(mapfile);
  transfer.WriteCorrespondMesh();

  //transfer
  transfer.Transfer(pose_body);
  transfer.WriteDeformedResult();

  return 0;
}
