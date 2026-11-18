#!/usr/bin/env bash

# Pinned commits of the external repositories that ./scripts/download-core.sh
# clones into ./distfiles/ for the non-Nix ./build.sh flow.
#
# Consumed by:
#   ./scripts/download-core.sh         (checks out exactly these commits)
#   ./scripts/update-core-versions.sh  (advances them to the upstream tips)
#
# Do not edit by hand; run
#   ./scripts/update-core-versions.sh            # show current vs upstream
#   ./scripts/update-core-versions.sh --update   # rewrite this file
#
# Format: one repository per line, fields separated by '|':
#
#   <name>|<clone url>|<branch>|<40-hex commit>|<commit date>
#
# <branch> is the upstream branch the commit was taken from; it is metadata
# used by ./scripts/update-core-versions.sh (which advances the pin to that
# branch tip) and a hint for humans. The commit is what is actually checked
# out, so a pin can be an older commit than the branch tip.
#
# <commit date> is the committer date (YYYY/MM/DD) of the pinned commit, like
# the dated pin comments in nixpkgs/gpt-lehner.nix. Informational only.
#
# The same sources are pinned for the Nix flow by nixpkgs/gpt-lehner.nix,
# nixpkgs/grid-lehner.nix and nixpkgs/cps.nix; advance those with
#   ./nixpkgs/upgrade.py --update gpt grid cps

CORE_REPOS=(
    # Grid family (built by ./scripts/grid*.sh).
    "Grid|https://github.com/paboyle/Grid.git|develop|3232b438f2cf729e515c185e2dbf29d7cc276bb7|2026/10/05"
    "Grid-clehner|https://github.com/lehner/Grid.git|feature/gpt|0d0774789f50eb3c453f23299de39eb81ceefce1|2026/10/25"
    "Grid-ljin|https://github.com/lehner/Grid.git|feature/gpt|0d0774789f50eb3c453f23299de39eb81ceefce1|2026/10/25"
    "Grid-tblum|https://github.com/tblum2/Grid.git|develop|f6917aeba3f56954329f6140cedcb24e4c51dfc1|2026/08/26"
    # GPT family (built by ./scripts/gpt.sh).
    "gpt|https://github.com/lehner/gpt.git|master|f0c024342e6440d62e6d258876a2e981a2c50072|2026/11/16"
    "gpt-ljin|https://github.com/lehner/gpt.git|master|f0c024342e6440d62e6d258876a2e981a2c50072|2026/11/16"
    # Hadrons family (built by ./scripts/hadrons*.sh).
    "Hadrons|https://github.com/aportelli/Hadrons.git|develop|1d64626aaa84052d33880bb028ef28de1d8f3047|2026/10/15"
    "Hadrons-tblum|https://github.com/tblum2/Hadrons.git|develop|9dc558bdfa3fa0bf970ded4a52b9053798444d39|2026/08/24"
    "Hadrons-vaishakhiN|https://github.com/vaishakhi123/Hadrons-tblum.git|develop|23cc53d58997eb3d43e85724f438d167a082ce7f|2026/08/19"
    # CPS (built by ./scripts/cps.sh; qlat-cps links libcps.a).
    "CPS|https://github.com/RBC-UKQCD/CPS_public.git|master|d3c8dd5e8a3ea6a315fd2fba963bf32585ed6331|2023/04/14"
)
