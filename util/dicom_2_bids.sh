#!/bin/bash

dicom_path="/home/pw0032/Data/FUS-RCT/dicom/sub-018-RCT/ses-2-00/"

sub_name=sub-018-RCT
session=ses-2-00

for d in $(find $dicom_path -type d); do
  if [ -z "$(find "$d" -mindepth 1 -type d)" ]; then
    dcm2niix -f "${sub_name}_${session}_%d_%s" -p y -z y -o $dicom_path $d
  fi
done


