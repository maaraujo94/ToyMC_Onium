#!/bin/sh

com="qqbarMC.C"

root -l <<EOF
gSystem->Load("/home/mariana/local/lib/libLHAPDF.so");
.include /home/mariana/local/include
.x $com
.q
EOF

echo "ending of file here: [state]_beta[val]"
read fname

mv MC_res.root MC_res_$fname.root
