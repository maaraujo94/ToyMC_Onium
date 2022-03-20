#!/bin/sh

com="ZMC.C"

root -l <<EOF
gSystem->Load("/home/mariana/local/lib/libLHAPDF.so");
.include /home/mariana/local/include
.x $com
.q
EOF
