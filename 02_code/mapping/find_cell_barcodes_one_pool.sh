DATA_DIR="$HOME/car_t_sc/01_data"
FASTQS_DIR="$DATA_DIR/raw/raw/2024-06-17_Maik_Luu_24054SC_raw_FASTQ/FASTQ_by_lib"
OUTPUT="$DATA_DIR/processed/count_cart_receptor"
SAM_INPUT="$DATA_DIR/processed/bowtie2_mapped_carts"
SCRIPT_DIR="$HOME/car_t_sc/02_code/mapping"

output_file=$OUTPUT/P1_GEX.csv
sam_file=$SAM_INPUT/P1_GEX.sam
jobname="Extract_Cellbarcode_P1_GEX"
file=$FASTQS_DIR/24054SC_Luu_P1_D0_GEX_B11_TT_R1.fq.gz

 sbatch --job-name="$jobname" \
           --error="${jobname}_err.txt" \
           --output="${jobname}_output.out" \
           --cpus-per-task=2 \
           --mem-per-cpu=2G \
           --wrap="python3 ${SCRIPT_DIR}/find_cell_barcodes.py -sam '$sam_file' -R1 "$file" -out '$output_file'"