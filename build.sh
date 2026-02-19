#!/bin/bash

mkdir -p $PREFIX/bin
cp gmaamp $PREFIX/bin/
chmod +x $PREFIX/bin/gmaamp

TARGET_DIR=$PREFIX/share/GMA
mkdir -p $TARGET_DIR

cp -r scripts/* $TARGET_DIR/
cp -r db/* $TARGET_DIR/