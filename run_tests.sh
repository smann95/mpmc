#! /bin/bash

newest_commit=$(git rev-parse HEAD)


for commit in $(git rev-list master)
do
  git checkout -b "${commit}" ${commit}
  bash compile.sh
  exit
done

#git checkout master
#git branch | grep -v "master" | grep -v "mpmc_tester" | xargs git branch -D
