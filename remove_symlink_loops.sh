#!/usr/bin/env bash

git_root="$(git rev-parse --show-toplevel)"

for dir in configs containers data scripts workflows; do
	rm -f "$git_root/$dir/$dir"
done
