#!/usr/bin/env bash

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


### For some obscure reason, nfcore/rnaseq decides to rename samples like MCF10A_T?_REPN* to MCF10A_REP?_T1*
### This script tries to correct file names and content (the latter only for text files)

set -e
set -o pipefail
set -u

if [ $# -lt 3 ]; then
  1>&2 echo "Usage: $0 root_dir query replacement [--replace]"
  1>&2 echo "Example: $0 data/output/nfcore_rnaseq MCF10A_REP1_T1 MCF10A_T1_REP1"
  exit 1
fi

root_dir="$1"
query="$2"
replacement="$3"

dry_run=true
if [[ $# -eq 4 && "$4" == '--replace' ]]; then
  dry_run=false
fi

function rename_file {
  src="$1"
  dest="${src//"$query"/"$replacement"}"

  1>&2 echo "mv $src $dest"

  if [ "$dry_run" == false ]; then
    mv "$src" "$dest"
  fi
}

function replace_text {
  fpath="$1"
  ftype="$(file -i "$fpath" | cut -d' ' -f2)"

  if [[ "$ftype" == 'text/plain;' ]]; then
    if perl -p -e '$M += s|$ENV{query}|$ENV{replacement}|; END{exit 1 unless $M>0}' "$fpath" > /dev/null; then
      if [ "$dry_run" == false ]; then
        1>&2 echo "Replacing pattern in $fpath"
        perl -pi -e 's|$ENV{query}|$ENV{replacement}|g' "$fpath"
      else
        1>&2 echo "One or more replacement would be made to file $fpath"
      fi
    fi
  fi
}

export -f rename_file replace_text
export query replacement dry_run

find "$root_dir" -type f -name "*$query*" -print0 |
  xargs -0 -I {} bash -c 'rename_file "$@"' _ {}

find "$root_dir" -type f -print0 |
  xargs -0 -I {} bash -c 'replace_text "$@"' _ {}
