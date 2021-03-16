//
// Created by jie zou on 2021/3/16.
//

#include <gtest/gtest.h>

#include <iostream>

class CommonVars: public ::testing::Test {
 protected:
  void SetUp() override {
    str_ = "gtesting...";
  }

 std::string str_;
};

TEST_F(CommonVars, testtest) {
  std::cout << str_ << std::endl;
}