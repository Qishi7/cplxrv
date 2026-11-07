# Publication to the owner's existing Qishi7/cplxrv repository

The owner explicitly requested replacing the contents of Qishi7/cplxrv rather
than creating a new fork. The original main commit is
7ef7f9e3decb6294e6774d08d9444768f6ad9f58 and remains in repository history.
The replacement installs package cplxshrink 0.1.0.9000; it does not retain the
old cplxrv public interface. The old package can be recovered from that commit.
No force-push, upstream write or upstream pull request is needed.

Replacement source comes from local development handoff commit
d004dff903f948a61b863fbf5f5e496f1a85a861, whose tested package code was
852bff6227a8b6c45f91625eba5dc8f0236d6c15. README installation instructions,
DESCRIPTION repository URLs and supervisor publication instructions are updated.
The four certified sampler/binding files and reference fixtures remain byte-identical
to frozen manuscript commit d985957c4254593ae2cc5c449d6dc3d108dd1b5b.
The new tree excludes the old tracked compiled binaries, .Rhistory, .Rproj.user
and .DS_Store files. Their previous history is not rewritten.

Historical baseline, checkpoint and final validation reports describe the initial
local run before GitHub upload. They remain valid historical evidence; their
statements that no fork/push had occurred refer to that earlier handoff.
This repository is not represented as a GitHub fork of upstream.

Publication status: SUCCESS. Replacement commit
e5bc100d9767a510a36a84392cdd3a5d9d9ac947 was pushed to main at
https://github.com/Qishi7/cplxrv. The remote main SHA was independently read
and matched this commit. The previous main commit remains its direct parent.
An initial Git transport attempt returned HTTP 400; retrying with a per-command
HTTP buffer/version setting succeeded without force-pushing.
GitHub CLI authorization was completed by the owner as Qishi7. No token or
credential is stored in this repository.

Fresh replacement-tree validation logs are stored under review/logs/publication/.
The cross-platform workflow is present; remote results require a successful push.

Replacement-tree validation completed before upload:
- 444 expectations, 0 failures/errors/warnings/skips: review/logs/publication/tests.log.
- R CMD build: succeeded; review/logs/publication/build.log.
- Standard R CMD check: Status OK, 0 errors/warnings/notes, including examples,
  both original test scripts and PDF manual: review/logs/publication/check.log.
- git diff against frozen manuscript commit for protected files/fixtures: exit 0;
  review/logs/publication/core-diff.log.
Only destination metadata/documentation and packaging ignore rules differ from
previously tested package code. GitHub CLI was downloaded from the official cli/cli
release and its published SHA-256 checksum verified; no system permissions changed.

## Remote CI status

The workflow was triggered by replacement commit e5bc100d9767a510a36a84392cdd3a5d9d9ac947:
https://github.com/Qishi7/cplxrv/actions/runs/42383210930
All four jobs failed before any runner steps started because of an account-level
GitHub Actions restriction. See the run annotations for details; the owner was
informed directly. No remote compilation, tests or platform checks ran, so Linux,
Windows and R 4.1 remain NOT VERIFIED. This is not reported as a package-test
failure or a CI pass. No account/billing settings were changed.
After resolving the account restriction, rerun the workflow to obtain cross-platform
validation. This evidence-only report update skips CI to avoid redundant jobs.

The complete installable source ZIP, source tarball, baseline diff, checksums and
local test/check reports are regenerated for this publication. Historical source
commit IDs identify the earlier local validation, while the replacement commit
identifies the published package. Code remains identical to the local tested
sampler/interfaces; only destination metadata and documentation were updated.
