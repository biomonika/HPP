set -e
set -x

#PAN011 #HG06804
PAN011_nanopore_list=(
"s3://human-pangenomics/submissions/37162c8a-a266-47d0-83e7-5b190e13bc69--UCSC_WashU_AFR/HG06804/nanopore/03_15_22_R941_HG06804_1.tar" 
"s3://human-pangenomics/submissions/37162c8a-a266-47d0-83e7-5b190e13bc69--UCSC_WashU_AFR/HG06804/nanopore/03_15_22_R941_HG06804_2.tar"
"s3://human-pangenomics/submissions/37162c8a-a266-47d0-83e7-5b190e13bc69--UCSC_WashU_AFR/HG06804/nanopore/03_15_22_R941_HG06804_3.tar" 
"s3://human-pangenomics/submissions/37162c8a-a266-47d0-83e7-5b190e13bc69--UCSC_WashU_AFR/HG06804/nanopore/03_15_22_R941_HG06804_4.tar" 
"s3://human-pangenomics/submissions/37162c8a-a266-47d0-83e7-5b190e13bc69--UCSC_WashU_AFR/HG06804/nanopore/03_15_22_R941_HG06804_5.tar" 
"s3://human-pangenomics/submissions/37162c8a-a266-47d0-83e7-5b190e13bc69--UCSC_WashU_AFR/HG06804/nanopore/03_15_22_R941_HG06804_6.tar"
"s3://human-pangenomics/submissions/37162c8a-a266-47d0-83e7-5b190e13bc69--UCSC_WashU_AFR/HG06804/nanopore/03_15_22_R941_HG06804_7.tar" 
)


for tar in ${PAN011_nanopore_list[@]}; do

  #echo $tar
  prefix="s3:\/"
  suffix=${tar//$prefix} #remove AWS prefix

  tar_name="$(basename -- $suffix)"
  echo $tar_name

  #check if tar file already exists in the right location and extracted
  tar_file_name_without_extension="${tar_name%.*}"
  extracted_raw_folder="/data/mcechova/guppy/${tar_file_name_without_extension}"
  mkdir -p ${extracted_raw_folder}

  download_folder="/public/groups/migalab/washu_pedigree/nanopore/raw"
  downloaded_tar=${download_folder}/${tar_name}

  if [ -f ${downloaded_tar} ] 
  then
      echo "Directory with the downloaded tar file already exists." 
      tar xvf ${downloaded_tar} -C ${extracted_raw_folder}
      rm ${downloaded_tar} #delete the tar after an extraction
  else
      echo "Directory with the extracted tar file DOES NOT exist."
      aws --no-sign-request s3 cp ${tar} ${tar_name}
      echo "Create a directory with .fast from the raw data tar"
      tar xvf ${tar_name} -C ${extracted_raw_folder}
      rm ${tar_name} #delete the tar after an extraction
  fi
  #basecall the given tar file
  ./basecallNanoporeFromTar_GPU4-7.sh ${extracted_raw_folder}
done

echo "Done."

