#!/bin/bash
# Splits full data table into mouse specific files
full_data_path=../data/sequence_data/m293_mouse_flu_bigtable_isotype

mkdir -p ../processed_data/mouse_specific_data_files/

# Comma separated string of unique mouse ids
#ids=$(grep -o '\s[0-9]*-[0-9]*\s' $full_data_path | sort | uniq | tr -d '\t' | tr '\n' ',')
#ids=${ids::-1}
ids=8-1,8-2,8-3,8-4,8-5,8-6,8-7,8-8,8-9,8-10,16-1,16-2,16-4,16-5,16-6,16-7,16-8,16-9,16-10,24-1,24-2,24-3,24-4,24-5,40-1,40-2,40-3,40-4,40-5,40-6,40-7,40-8,40-9,40-10,56-1,56-2,56-3,56-4,56-5,56-6,56-7,56-8,56-9,56-10

for mouse_id in $(echo $ids | sed "s/,/ /g")
do
    echo "Exporting sequences of mouse $mouse_id"
    output_file_path=../processed_data/mouse_specific_data_files/$mouse_id.csv
    
    # Comma separated string of unique mouse ids
    head -1 $full_data_path | tr "\t" "," > $output_file_path

    # Extract lines for individual mouse, replacing tabs with comma
    grep "\s${mouse_id}\s" $full_data_path | tr "\t" "," >> $output_file_path
    
done

echo "Done."