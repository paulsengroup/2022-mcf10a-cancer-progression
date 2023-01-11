#!/usr/bin/env bash

git_root="$(git rev-parse --show-toplevel)"

for dir in bin bin configs containers data workflows; do
	rm -f "$git_root/$dir/$dir"
done
