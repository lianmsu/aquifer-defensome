#!/bin/bash
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

minced -spacers ${HOME_ROOT}/03GW/12minced/1cas_seq/cas_seq.fna bin.crisprs bin.gff
