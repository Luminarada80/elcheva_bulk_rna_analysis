MAX_JOBS_IN_QUEUE=5

BASE_DIR="/gpfs/Labs/Uzun/SCRIPTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/elcheva_bulk_rna_analysis"

SAMPLE_NAMES=(
    # "K562_shCNTL_Rep1"
    "K562_shCNTL_Rep2"
    "K562_shCNTL_Rep3"

    "K562_shNYN3_Rep1"
    "K562_shNYN3_Rep2"
    "K562_shNYN3_Rep3"

    "K562_shNYN4_Rep1"
    "K562_shNYN4_Rep2"
    "K562_shNYN4_Rep3"

    "KG1A_shCNTL_Rep1_S18"
    "KG1A_shCNTL_Rep2_S20"
    "KG1A_shCNTL_Rep3_S22"

    "KG1A_shNYN4_Rep1_S19"
    "KG1A_shNYN4_Rep2_S21"
    "KG1A_shNYN4_Rep3_S23"

    "SKNAS_shCNTL_Rep1"
    "SKNAS_shCNTL_Rep2"
    "SKNAS_shCNTL_Rep3"

    "SKNAS_shNYN4_Rep1"
    "SKNAS_shNYN4_Rep2"
    "SKNAS_shNYN4_Rep3"
)

for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do

    LOGDIR=${BASE_DIR}/LOGS/${SAMPLE_NAME}
    mkdir -p $LOGDIR

    sbatch \
        --export=ALL,SAMPLE_NAME="$SAMPLE_NAME" \
        --job-name="RNAseqpipeline_${SAMPLE_NAME}" \
        --output="${LOGDIR}/${SAMPLE_NAME}.log" \
        --error="${LOGDIR}/${SAMPLE_NAME}.err" \
        "${BASE_DIR}/RNAseqpipeline.sh"

done