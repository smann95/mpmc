#! /bin/bash
#set -e so if a command fails we immediately exit

set -e
newest_commit=$(git rev-parse HEAD)

mkdir -p compile_logs
mkdir -p run_logs
test_input_dir="tests"

for commit in $(git rev-list master | head)
do
  rm -rf build
  git checkout -b "${commit}" ${commit}
  bash compile.sh &> compile_logs/compile_${commit}.log
  build/mpmc bs.inp
  exit
done

git checkout test_suite
git branch | grep -v "master" | grep -v "test_suite" | xargs git branch -D
