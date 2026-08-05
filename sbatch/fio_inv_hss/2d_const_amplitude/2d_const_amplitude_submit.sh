#!/bin/bash
# =============================================================================
# Sequential SLURM job submission script
# Runs 2d_const_amplitude_cg.sh first, then submits 2d_const_amplitude_indep.sh
# only if the first job exits successfully (exit code 0).
# =============================================================================

set -euo pipefail

SCRIPT1="2d_const_amplitude_cg.sh"
SCRIPT2="2d_const_amplitude_indep.sh"

# ---- Check that the script files exist ----
if [[ ! -f "$SCRIPT1" ]]; then
    echo "Error: First script not found: $SCRIPT1"
    exit 1
fi

if [[ ! -f "$SCRIPT2" ]]; then
    echo "Error: Second script not found: $SCRIPT2"
    exit 1
fi

# ---- Submit the first job ----
echo "Submitting first job: $SCRIPT1 ..."
JOB1_OUTPUT=$(sbatch "$SCRIPT1" 2>&1)

# Extract Job ID (sbatch output format: Submitted batch job 12345)
JOB1_ID=$(echo "$JOB1_OUTPUT" | awk '/Submitted batch job/{print $4}')

if [[ -z "$JOB1_ID" ]]; then
    echo "Error: Failed to submit the first job"
    echo "Output: $JOB1_OUTPUT"
    exit 1
fi

echo "  -> Job ID: $JOB1_ID"

# ---- Submit the second job, dependent on the first job succeeding ----
echo "Submitting second job: $SCRIPT2 (depends on Job $JOB1_ID succeeding)..."
JOB2_OUTPUT=$(sbatch --dependency=afterok:"$JOB1_ID" "$SCRIPT2" 2>&1)

JOB2_ID=$(echo "$JOB2_OUTPUT" | awk '/Submitted batch job/{print $4}')

if [[ -z "$JOB2_ID" ]]; then
    echo "Error: Failed to submit the second job"
    echo "Output: $JOB2_OUTPUT"
    exit 1
fi

echo "  -> Job ID: $JOB2_ID"

# ---- Summary ----
echo ""
echo "========================================"
echo "Submission complete:"
echo "  [1] $SCRIPT1  -> Job $JOB1_ID"
echo "  [2] $SCRIPT2  -> Job $JOB2_ID (afterok:$JOB1_ID)"
echo "========================================"
echo ""
echo "Check queue:  squeue -u \$USER"
echo "Check deps:   squeue -u \$USER -o '%.10i %.20j %.2t %.10M %.20R %E'"