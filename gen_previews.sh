#!/bin/sh

SCENE_DIR=assets/scenes
PREVIEW_DIR=assets/previews

./raytrace $SCENE_DIR/test1.scene $PREVIEW_DIR/test1.jpg &&
./raytrace $SCENE_DIR/test2.scene $PREVIEW_DIR/test2.jpg &&
./raytrace $SCENE_DIR/spheres.scene $PREVIEW_DIR/spheres.jpg &&
./raytrace $SCENE_DIR/table.scene $PREVIEW_DIR/table.jpg &&
./raytrace $SCENE_DIR/SIGGRAPH.scene $PREVIEW_DIR/SIGGRAPH.jpg 

