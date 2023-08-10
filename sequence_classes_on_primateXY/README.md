##This set of script annotates sequence classes on the X and Y chromosome of primates species

The scripts should be run in the following order:
1) run prepare_ampliconic_regions.sh (the script analyzes the results of the intrachromosomal similarity as revealed by blast, and filters them) )
2) run create_sequence_classes_Y.sh (the script generates both bed files for further downstream analysis, as well as files suitable for plotting with circos)
3) run create_sequence_classes_X.sh (the script generates both bed files for further downstream analysis, as well as files suitable for plotting with circos)
4) manually adjust XTR for both X and Y annotations
