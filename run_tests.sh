#! /bin/bash

newest_commit=$(git rev-parse HEAD)


for commit in $(git rev-list master | head)
do
  git checkout -b "${commit}" ${commit}
  bash compile.sh
done

git checkout test_suite
git branch | grep -v "master" | grep -v "test_suite" | xargs git branch -D
