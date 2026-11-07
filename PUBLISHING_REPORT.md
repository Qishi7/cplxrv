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

Upload status at preparation: NOT YET PUSHED. The connector returned HTTP 403
"Resource not accessible by integration" for Git tree writes; command-line Git
had no credentials. An official GitHub CLI login is being prepared.
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
