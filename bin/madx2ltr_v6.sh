#!/bin/bash
if [ $# != 1 ]; then
echo "Script makes lifetrac task file from madx output"
echo "must have out.lattice, esave, out.strong in current dir."
echo "by A.Valishev (valishev@fnal.gov)"
echo " "
echo "usage: $0 task_file"
exit 0;
fi

perl ~/lifetrac/bin/madx2ltr_v6.pl out.lattice esave out.strong $1
