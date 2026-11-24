#!/usr/bin/env bash
# Run a Snakemake target inside whole-node SLURM jobs, one after another.
#
# Each job drives Snakemake with profiles/local on the GPUs it was given. The
# job submits its successor first, held until this job ends. Before the time
# limit, SLURM sends SIGINT to this script, which passes it to Snakemake and
# its jobs. Snakemake exits without a stale lock, and the successor reruns the
# jobs that were cut off. When Snakemake ends on its own, the job cancels its
# successor and the chain stops.
#
# usage: scripts/node_chain.sh <sbatch options> -- <snakemake arguments>
#   e.g. scripts/node_chain.sh -p gpuh200 --gres=gpu:8 -c 64 -t 1-0 -- \
#          train --configfile config/<run_id>/config.yaml
# Request the GPUs explicitly (--gres=gpu:N). Where cgroups hide GPUs a job did
# not request, an exclusive node alone leaves nvidia-smi with none to list.
# CHAIN_GRACE (seconds, default 300) is how long before the limit SIGINT comes.
# SLURM may send it up to a minute late.
set -uo pipefail

if [[ -z ${SLURM_JOB_ID:-} ]]; then
	sbatch_options=()
	while [[ $# -gt 0 && $1 != -- ]]; do
		sbatch_options+=("$1")
		shift
	done
	[[ ${1:-} == -- ]] || { echo "usage: $0 <sbatch options> -- <snakemake arguments>" >&2; exit 2; }
	shift
	mkdir -p logs/chain
	CHAIN_SCRIPT=$(realpath "$0")
	CHAIN_SBATCH_OPTIONS=$(printf '%s\n' "${sbatch_options[@]}")
	export CHAIN_SCRIPT CHAIN_SBATCH_OPTIONS
	exec sbatch "${sbatch_options[@]}" \
		--job-name="${CHAIN_NAME:-node-chain}" \
		--output=logs/chain/%x.%j.out \
		--signal="B:INT@${CHAIN_GRACE:-300}" \
		"$CHAIN_SCRIPT" "$@"
fi

cd "$SLURM_SUBMIT_DIR" || exit 1
mapfile -t sbatch_options <<<"$CHAIN_SBATCH_OPTIONS"

successor=$(sbatch --parsable "${sbatch_options[@]}" \
	--job-name="$SLURM_JOB_NAME" \
	--output=logs/chain/%x.%j.out \
	--signal="B:INT@${CHAIN_GRACE:-300}" \
	--dependency="afterany:$SLURM_JOB_ID" \
	"$CHAIN_SCRIPT" "$@")
echo "[chain] job $SLURM_JOB_ID on $(hostname), successor $successor"

# `pixi run` does not pass SIGINT on, so Snakemake runs from the activated env.
eval "$(pixi shell-hook)"

n_gpus=$(nvidia-smi -L | wc -l)
# Job control gives Snakemake its own process group. Without it, bash also
# starts a background command with SIGINT ignored.
set -m
snakemake "$@" \
	--profile profiles/local \
	--resources "gpu=$n_gpus" \
	--cores "${SLURM_CPUS_PER_TASK:-$(nproc)}" &
driver=$!
set +m

interrupted=false
trap 'interrupted=true; echo "[chain] time limit near, stopping Snakemake"; kill -INT -- "-$driver"' INT
while kill -0 "$driver" 2>/dev/null; do
	wait "$driver"
	status=$?
done
echo "[chain] Snakemake exited with $status, interrupted=$interrupted"

if ! $interrupted; then
	scancel "$successor"
	echo "[chain] chain ends, cancelled successor $successor"
fi
exit "$status"
