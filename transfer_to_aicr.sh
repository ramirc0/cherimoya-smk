#!/usr/bin/env bash
set -euo pipefail

# One-way rclone copies between this workdir and AICR. `copy` never deletes.
#   ./transfer_to_aicr.sh push [rclone flags]   workflow, config, refs, staged inputs, notes -> AICR
#   ./transfer_to_aicr.sh pull [rclone flags]   results, logs, benchmarks -> here
# Pass --dry-run to preview.

cd "$(dirname "$0")"
cmd="${1:?usage: transfer_to_aicr.sh <push|pull> [rclone flags]}"
shift

remote="AICR:/work/umass/genomeatlas_umassmed/cherimoya-smk"
opts=(--transfers 16 --checkers 32 --sftp-concurrency 64 --multi-thread-streams 4 -P)

case "$cmd" in
  push)
    # Envs are rebuilt on AICR. --copy-links sends the files behind the
    # resources/ symlinks (refs, motifs, the staged inputs), not the links.
    # rclone applies the first matching filter, so the env excludes precede
    # the caller's filters.
    rclone copy . "$remote" "${opts[@]}" --copy-links \
      --filter '- .pixi/**' \
      --filter '- __pycache__/**' \
      "$@" \
      --filter '+ /workflow/**' \
      --filter '+ /config/**' \
      --filter '+ /profiles/**' \
      --filter '+ /tests/**' \
      --filter '+ /docs/**' \
      --filter '+ /resources/refs/**' \
      --filter '+ /resources/folds/**' \
      --filter '+ /resources/motifs/**' \
      --filter '+ /resources/chip-nexus/**' \
      --filter '+ /resources/procap/**' \
      --filter '+ /resources/chip-exo/**' \
      --filter '+ /resources/*.bed' \
      --filter '+ /pixi.toml' \
      --filter '+ /pixi.lock' \
      --filter '+ /pyproject.toml' \
      --filter '+ /pytest.ini' \
      --filter '+ /README.md' \
      --filter '+ /.claude/CLAUDE.md' \
      --filter '+ /.claude/handoffs/**' \
      --filter '+ /.gitignore' \
      --filter '+ /transfer_to_aicr.sh' \
      --filter '- **'
    ;;
  pull)
    rclone copy "$remote" . "${opts[@]}" "$@" \
      --filter '+ /results/**' \
      --filter '+ /logs/**' \
      --filter '+ /benchmarks/**' \
      --filter '- **'
    ;;
  *)
    echo "unknown command: $cmd (use push | pull)" >&2
    exit 1
    ;;
esac
