#!/bin/sh

com="gammaMC_pT.C"

root -l <<EOF
gSystem->Load("/home/mariana/local/lib/libLHAPDF.so");
gSystem->AddIncludePath("/home/mariana/local/include")
.x $com
.q
EOF
