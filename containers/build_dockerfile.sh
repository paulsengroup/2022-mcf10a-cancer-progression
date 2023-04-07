# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

set -e
set -o pipefail
set -u

if [ $# -ne 1 ]; then
  1>&2 echo "Usage: $0 Dockerfile"
  1>&2 echo "Example: $0 utils__v1.0.1.Dockerfile"
  exit 1
fi

git_root="$(git rev-parse --show-toplevel)"

dockerfile="$1"
dockerfile_name="$(basename "${dockerfile[*]}" .Dockerfile)"
name="$(echo "$dockerfile_name" | sed -E 's/(.*)__.*/\1/')"
version="$(echo "$dockerfile_name" | sed -E 's/.*__v(.*)/\1/')"

sudo docker build \
  -t "$name:latest" \
  -t "$name:$version" \
  -f "$dockerfile" \
  --build-arg="CONTAINER_VERSION=$version" \
   "$git_root"

if ! command -v apptainer &> /dev/null; then
  2>&1 echo 'Unable to find apptainer in your PATH'
  2>&1 echo 'Skipping apptainer build step!'
  exit 0
fi

sif="$git_root/containers/cache/ghcr.io-paulsengroup-2022-mcf10a-cancer-progression-$name-$version.img"
sudo apptainer build -F "$sif" "docker-daemon://$name:$version"

uid="$(id -u)"
gid="$(id -g)"

sudo chown "$uid:$gid" "$sif"

set -x
apptainer exec "$sif" sh -c 'echo Hi!' > /dev/null
