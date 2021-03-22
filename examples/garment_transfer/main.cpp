//
// Created by jie zou on 2020/8/7.
//
#include "garment_transfer.h"

int main() {
  const std::string cloth     = DATA_PATH"/cloth.obj";
  const std::string body      = DATA_PATH"/dress_body.obj";
  const std::string mapfile   = DATA_PATH"/cloth2body.map";

  const std::string pose_body = DATA_PATH"/pose_body.obj";

  const std::string cor_cloth = DATA_PATH"/cor_colth.obj";
  const std::string de_cloth  = DATA_PATH"/de_colth.obj";
  const std::string fit_cloth = DATA_PATH"/fit_colth.obj";

  pmp_jie::GarmentTransfer transfer;

  //correspond
  transfer.Load(cloth);
  transfer.LoadTarget(body);
  transfer.Deform(mapfile);
  transfer.WriteCorrespondMesh(cor_cloth);

  //transfer
//  transfer.Transfer(pose_body);
  transfer.PositionTransfer(pose_body, {2344});
  transfer.WriteDeformedResult(de_cloth);
  transfer.Fit(de_cloth);
  transfer.WriteGarment(fit_cloth);

  return 0;
}
