#!/usr/bin/env bash
set -euo pipefail

# Repair fMRIPrep batch outputs created through sudo Docker:
# - fix ownership and permissions on the full tree
# - replace FreeSurfer symlinks inside sourcedata/freesurfer with real files

ROOT_DIR="${1:-/home/sskgroup/Documents/fmriprep_ch}"
TARGET_USER="${SUDO_USER:-sskgroup}"
TARGET_GROUP="$(id -gn "$TARGET_USER" 2>/dev/null || echo sskgroup)"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

require_root() {
  if [[ ${EUID:-$(id -u)} -ne 0 ]]; then
    log "Re-running with sudo so ownership can be fixed..."
    exec sudo -E bash "$0" "$@"
  fi
}

flatten_symlinks_in_tree() {
  local fs_root="$1"
  local link target
  local changed pass=1

  log "Flattening symlinks under: $fs_root"

  while :; do
    changed=0

    while IFS= read -r -d '' link; do
      target="$(readlink -f -- "$link" || true)"

      if [[ -z "$target" || ! -e "$target" ]]; then
        log "WARN: broken link removed: $link"
        rm -f -- "$link"
        changed=1
        continue
      fi

      # Replace the symlink itself with a real file or directory copy.
      rm -f -- "$link"
      cp -a -- "$target" "$link"
      changed=1
    done < <(find "$fs_root" -type l -print0)

    [[ $changed -eq 0 ]] && break
    log "Pass $pass complete for $fs_root"
    pass=$((pass + 1))
  done
}

main() {
  require_root "$@"

  if [[ ! -d "$ROOT_DIR" ]]; then
    log "ERROR: root directory does not exist: $ROOT_DIR"
    exit 1
  fi

  log "Using root directory: $ROOT_DIR"
  log "Setting ownership to: $TARGET_USER:$TARGET_GROUP"
  chown -R "$TARGET_USER:$TARGET_GROUP" "$ROOT_DIR"

  log "Setting permissions to user-writable and group/world-readable"
  chmod -R u+rwX,go+rX "$ROOT_DIR"

  while IFS= read -r -d '' fs_root; do
    flatten_symlinks_in_tree "$fs_root"
  done < <(find "$ROOT_DIR" -type d -path '*/sourcedata/freesurfer' -print0)

  log "Re-applying ownership/permissions after link flattening"
  chown -R "$TARGET_USER:$TARGET_GROUP" "$ROOT_DIR"
  chmod -R u+rwX,go+rX "$ROOT_DIR"

  log "Done. Remaining symlinks under root: $(find "$ROOT_DIR" -type l | wc -l)"
}

main "$@"
