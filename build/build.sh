#!/bin/sh
set -e

mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/db"
cp -r bin/* "${PREFIX}/bin/"
cp -r db/* "${PREFIX}/db/"

