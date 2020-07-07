#!/bin/sh

./raytrace test1.scene test1.jpg &&
./raytrace test2.scene test2.jpg &&
./raytrace spheres.scene spheres.jpg &&
./raytrace table.scene table.jpg &&
./raytrace SIGGRAPH.scene SIGGRAPH.jpg 

