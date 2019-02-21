#!/bin/bash
# A new directory is created named Phoenix-<commit>. Then, for config file 
# given as command arguments, a new project will be created in the test 
# directory with the config file added. Finally, the pipeline will be launched
# as a slurm batch job.
set -ue 

REPO=~/jetstream_pipelines/phoenix/

if (cd ${REPO} && git diff --exit-code); then
    echo "${REPO} repo looks clean..."
else
    echo 'ERROR! Uncommitted changes found in repo!' && exit 1
fi

COMMIT="$(cd ${REPO} && git rev-parse --short HEAD)"
DIRNAME="Phoenix-${COMMIT}"
mkdir -p "${DIRNAME}"

for CONFIG in $@; do
  FILENAME="$(basename $CONFIG)"
  PROJECT="${FILENAME%.*}"
  jetstream init "${DIRNAME}/${PROJECT}" -- --config "${CONFIG}" 
  sbatch -o "${DIRNAME}/%x-slurm-%j.out" -J "${PROJECT}" --wrap "exec jetstream pipelines phoenix --project ${DIRNAME}/${PROJECT}"
done

