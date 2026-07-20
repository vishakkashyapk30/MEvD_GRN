#!/bin/bash
# =============================================================================
# MEvD-GRN full training on Ada (research account).
# Maintains the required 1 GPU : 10 CPU ratio. The model is tiny and trains in
# well under an hour; the time budget is generous for preprocessing + reruns.
#   submit:  sbatch slurm/train.sh
#   watch :  tail -f logs/train_<jobid>.out
# =============================================================================
#SBATCH --job-name=MEvD-GRN
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=long
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/train_%j.out
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs
cd "$SLURM_SUBMIT_DIR"

# --- environment ----------------------------------------------------------
# Custom envs must live in $HOME (Ada software policy). Create once:
#   conda create -y -n mevd-grn python=3.10 && conda activate mevd-grn
#   pip install torch==2.1.0 --index-url https://download.pytorch.org/whl/cu118
#   pip install torch_geometric torch_scatter torch_sparse \
#       -f https://data.pyg.org/whl/torch-2.1.0+cu118.html
#   pip install -r requirements.txt && pip install -e .
source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate mevd-grn
module load u18/cuda/11.6 2>/dev/null || true

python -c "import torch; print('CUDA:', torch.cuda.is_available(), torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU')"

# --- preprocess (skip if already done) -----------------------------------
for ct in K562 ESC; do
  cfg="configs/$(echo $ct | tr 'A-Z' 'a-z').yaml"
  if [ ! -f "data/processed/$ct/rna_features_aligned.npy" ]; then
    echo "### Preprocessing $ct ###"
    python scripts/02_preprocess.py --config "$cfg"
  fi
done

# --- train both cell types -----------------------------------------------
echo "### Training K562 ###"
python scripts/03_train.py --config configs/k562.yaml --device cuda:0
echo "### Training ESC ###"
python scripts/03_train.py --config configs/esc.yaml --device cuda:0

# --- cross-cell-type transfer (K562 -> ESC) ------------------------------
echo "### Transfer eval ###"
python scripts/04_evaluate.py --config configs/k562.yaml \
  --checkpoint results/checkpoints/K562/final_model.pt \
  --transfer_config configs/esc.yaml --device cuda:0

echo "### DONE ###"
