// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Stores metadata like git commit of current build.

#pragma once

#include <iostream>
#include <sstream>

class Version {
 public:
  static std::string GetGitCommit() { return git_hash_; }
  static std::string GetGitBranch() { return git_branch_; }
  static std::string GetGitTags() { return git_tags_; }

  static inline const std::string git_hash_ = GIT_HASH;
  static inline const std::string git_branch_ = GIT_BRANCH;
  static inline const std::string git_tags_ = GIT_TAGS;
};