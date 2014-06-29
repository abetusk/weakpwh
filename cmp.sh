#!/bin/bash
#

g++ -g clipper.cpp sweep/advancing_front.cc sweep/cdt.cc sweep/sweep.cc sweep/sweep_context.cc common/shapes.cc weakpwh.cpp -o weakpwh
#g++ -g sweep/advancing_front.cc sweep/cdt.cc sweep/sweep.cc sweep/sweep_context.cc common/shapes.cc weakpwh.cpp -o weakpwh

