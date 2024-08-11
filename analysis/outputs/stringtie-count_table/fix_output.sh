source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg

STRINGTIE_OUTPUT_DIR=$out/stringtie-count_table


# ## rename .out to .gtf
# find $STRINGTIE_OUTPUT_DIR -type f -name "*.out" | while read file; do
#   mv "$file" "${file%.out}.gtf"
# done

# # rename {treatment}.{ID} to {treatment}_{ID}
# find $STRINGTIE_OUTPUT_DIR -type f -name "*.gtf" | while read file; do
#   mv "$file" "${base/.SRR/_SRR}.count-table.gtf"
# done

#mv all .gtf files to a collective folder, but keep individuals for .ctab info etc.
# mkdir -p $STRINGTIE_OUTPUT_DIR/gtf
# find $STRINGTIE_OUTPUT_DIR -type f -name "*.gtf" | while read file; do
#   mv $file $STRINGTIE_OUTPUT_DIR/gtf/$(basename $file)
# done

find $STRINGTIE_OUTPUT_DIR -type f -name "*.gtf" | while read file; do
  # mkdir -p $STRINGTIE_OUTPUT_DIR/gtf
  # mv $file $STRINGTIE_OUTPUT_DIR/gtf/$(basename $file)
  # echo $file
  # echo $(basename $file)
  output=$out/stringtie-prepDE/$(basename ${file%.count-table.gtf})
  mkdir -p $output
  cp $file $output/$(basename $file)
done