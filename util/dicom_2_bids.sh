#!/bin/bash

dicom_path="/mnt/msi/Data/naviFUS/dicom"
out_path="/mnt/msi/Data/naviFUS/nifti"

sub_name=sub-009-NAVI
session=ses-00

for d in $(find $dicom_path -type d); do
  if [ -z "$(find "$d" -mindepth 1 -type d)" ]; then
    dcm2niix -f "${sub_name}_${session}_%d_%s" -p y -z y -o $out_path/$sub_name/$session $d
  fi
done


