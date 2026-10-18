# Implementation state and planning

Updated 2026-10-18. Source baseline: `7a6a6b3` on `ben-bam-support`.
Read [architecture.md](architecture.md) for the current implementation and
[README.md](README.md) for context-maintenance rules.

## Current state

- Local Git history and `origin` were restored, and `ben-bam-support` was committed
  and pushed. Its published source commit is
  `7a6a6b397b5e8d3a0c18575dbe65a7e45a7bbed6`.
- The branch deliberately starts from upstream
  `687b9d3a385c11366fda54eef6e69df5f02ed358`, dated 2023-06-29.
  It has not incorporated later `master` changes.
- BAM parsing, test fixtures, local parser changes, operational notes, a Perl
  utility, and Docker changes are preserved as they existed before restoration.
  Source code has not been corrected as part of project-context setup.
- This setup adds `AGENTS.md`, `prompts/README.md`, `prompts/architecture.md`, and
  this plan. No previous project-local instruction or prompt files existed.
  Global/private startup configuration remains outside the repository scaffold.
- The four context files were reviewed for publication on 2026-10-18 under the
  user-authorized repository cleanup. The documentation publication commit
  follows the source preservation commit above; application source is unchanged.
- No human-authored development prompt or deep-record archive exists. There is
  no active implementation assignment beyond this documentation setup.

## Project-context setup plan

| Step | Status | Verification |
| --- | --- | --- |
| Read the development workflow and project-setup reference; inspect existing conventions. | Complete | Working tree initially clean; no existing scaffold to overwrite. |
| Review the full source, tests, configurations, packaging, and support scripts. | Complete | All 11 Python files and the Perl script read; fixture and archive inventories inspected. |
| Document architecture, interfaces, observed limitations, and Git restoration. | Complete | Four standard context files created with current facts separated from proposed work. |
| Verify routing, local links, factual checks, and scope of changes. | Complete | All four files re-read, 29 local links checked, both context indexes verified, and tracked source/fixtures unchanged. |

## How this directory was attached to origin

### Starting point and baseline decision

The directory originated from
`git@github.com:methylgrammarlab/epiread-tools.git`, but its `.git` directory had
been deleted and files had been edited outside Git. Fetching could recover
published history, not any unpushed commits lost with the deleted metadata.

The remote's default branch was `master`, at
`d97df3bfff3c7af8939f815214405f58a1b6bc3c` (2024-11-26).
The local tree lacked later parser changes and `requirements.txt`. Attaching
the snapshot directly to current `master` would have recorded those differences
as reversions/deletion. A comparison of upstream history found the closest
source-text match at `687b9d3`; this is an adopted restoration baseline, **not
proof of the original checkout commit**.

The user explicitly chose preserving the current files on a branch from that
historical commit, leaving integration with newer upstream changes for separate
work. `requirements.txt` is absent from the chosen base too, so the resulting
preservation commit contains no deletion of it.

### Backup and preparation

Before changing `.gitignore` or restoring Git metadata, the entire directory was
copied to this sibling directory on the original Mac:

```text
/Users/benbfly/Documents/epiread_tools.backup-20260916-094802
/Users/benbfly/Documents/epiread_tools.backup-20260916-094802.manifest.json
```

The backup contained 55 files totaling 204,451 bytes. Per-file SHA-256 values
verified the copy. The JSON manifest is outside the repository; it is recovery
evidence, not a file expected in other checkouts. Keep the backup until the user
authorizes its removal.

The local `.gitignore` originally contained only `deconvolution_models`.
Remote `master` had no root `.gitignore`, only `.idea/.gitignore`. The root file
was expanded to ignore Python bytecode, build/package outputs, pytest/coverage
outputs, local virtual environments, and `.DS_Store`. The optional blanket
`.idea/` exclusion was not added. Already tracked `dist/` packages and IDE files
were preserved; ignore rules do not untrack existing files.

Remote history was inspected in an isolated temporary bare clone at
`/tmp/epiread-origin-check.7sImi5`, initially shallow and then unshallowed.
That temporary inspection directory is not required for continued development.
The proposed remote branch name was confirmed absent before restoration.

### Git operations performed

These commands were run from the existing working directory after baseline
selection. They are a historical record, **not commands to rerun on this repo**.

```bash
git init -b ben-bam-support
git remote add origin git@github.com:methylgrammarlab/epiread-tools.git
git fetch origin
git reset --mixed 687b9d3a385c11366fda54eef6e69df5f02ed358
```

The mixed reset attached the new branch and index to the existing upstream
commit without replacing working files. No hard reset, checkout of upstream
files, merge, rebase, or force push was performed. This was ordinary Git history
attachment, not a `.git/info/grafts` or replace-ref operation.

After reviewing the differences and verifying file preservation:

```bash
git add -A
git commit -m 'Preserve local BAM support work and ignore generated files' \
  -m 'Restore local changes on upstream commit 687b9d3, including the ModbamReader, BAM fixtures, parsing changes, scripts, and Docker adjustments. Add ignore rules for Python caches, packaging outputs, test outputs, local environments, and macOS metadata.'
git push -u origin ben-bam-support
```

The commit was authored using the existing Git identity, Ben Berman. It contains
22 changed paths: five modifications and 17 additions, with 552 inserted and 17
removed text lines plus binary fixtures. The modified existing files were
`Dockerfile`, `em_utils.py`, `epiformat.py`, `epiparser.py`, and
`epireadToBedgraph.py`. Additions include `.gitignore`, `README.txt`, the Perl
script, `ModbamReader`, BAM-related fixtures/variants, JSON examples, and the BAM test.

### Verification and outcome

- Every original file except the authorized `.gitignore` edit matched the backup
  byte for byte, both after attachment and after the push.
- Every staged blob was independently compared with the local file bytes.
  There were no deletions or leftover nonignored untracked files.
- The preservation commit's parent is exactly `687b9d3`.
- Remote `ben-bam-support` was verified at `7a6a6b3`; local tracking was configured
  to `origin/ben-bam-support`, and the working tree was clean at completion.
- Remote `master` remained at `d97df3b`; no pull request or upstream merge was made.
- At setup inspection, the merge base with fetched `origin/master` was `687b9d3`:
  one local source commit versus eight later upstream commits. This count is a
  dated observation, not a permanent branch invariant.
- Application tests were not run during Git preservation. No further Git
  obstacle occurred after the user selected the historical baseline.

Branch: <https://github.com/methylgrammarlab/epiread-tools/tree/ben-bam-support>

## Validation baseline for project setup

Observed on this Mac on 2026-09-16:

| Check | Result and limit |
| --- | --- |
| Python syntax | All 11 Python source/test/setup files parsed with `ast.parse`; no imports or bytecode generation required. |
| Perl syntax | `perl -c ben_scripts/join_lambdas_thetas.pl` passed; no functional join test performed. |
| JSON examples | Both `tests/*.json` files parsed. `epiread_tools/sample_config.json` failed due to its trailing comma; its input paths also do not exist here. |
| Compressed fixtures | All seven `.gz` text fixtures decompressed successfully; corresponding `.tbi` files exist. This does not establish correct Tabix queries. |
| BAM fixtures | Both BAMs decompressed with BAM headers, `chr16` reference length 90,338,345, and two complete alignment blocks each. Both index files have BAI headers. This does not validate index contents or modification decoding. |
| Historical distributions | All seven package archives opened; none lists `ModbamReader`. They should not be used to test this branch's BAM code. |
| Runtime availability | `/usr/local/bin/python3` is Python 3.13.7; NumPy, pandas, SciPy, pysam, pytest, and Click were unavailable to that interpreter. No dependency installation was attempted. |
| Functional suite | Not run. Six existing tests identified by source inspection; no passing runtime baseline claimed. |
| Container | Not built. Required `deconvolution_models/` build-context directory is absent. |
| Documentation | Four files re-read; all 29 local links resolve; code fences and whitespace checked; both context indexes route all three prompt files. Exactly four new context files, with no staged or tracked-file changes. |

## Documentation publication review — 2026-10-18

- Reviewed all four context additions and the existing ignore rules. No extra
  ignore patterns were needed; authored documentation is retained.
- Confirmed one worktree, no nested Git repository, and the configured public
  `methylgrammarlab/epiread-tools` remote. Before publication, local and remote
  `ben-bam-support` both pointed to source preservation commit `7a6a6b3`.
- All 29 relative documentation links resolve; Markdown fences and whitespace
  checks pass. This is documentation-only work, so no application tests were run
  and the earlier functional-validation limitations remain applicable.
- Publish on the existing branch with an ordinary push, then verify the remote
  hash and clean working tree. No upstream integration is included.

## Proposed follow-up work

These are observed defects or validation gaps, not changes made by this setup
and not authorization to implement them. Confirm the next work item with the user.

| Priority / area | Observation | Focused verification for future work |
| --- | --- | --- |
| High: runnable baseline | Runtime dependencies are absent; package dependencies have no version bounds or tested Python range. | Choose an isolated environment, record versions, and run the six existing tests from the repo root before changes. |
| High: CLI/configuration | JSON list paths are passed to `.split()`; `intervals` does not match the declared `genomic_intervals` option; default format can become `None`; Click excludes `modbam`; documentation contradicts actual JSON precedence. | CLI tests for JSON-only, flags-only, overrides, all four formats, and output values; then align docs and examples. |
| High: BAM read semantics | A matrix row is allocated per call rather than read; relative index zero is discarded; unknown quality can become unmethylated; modification codes are not selected. | Tiny fixtures with multiple CpGs per read, first mapped CpG, thresholds/unknown quality, missing tags, and mixed modifications. Check row identity and exact states, not just coverage totals. |
| High: coordinates and sparse edge cases | Fetch start subtraction, inclusive end lookup, strict atlas membership, and implicit assumptions about CpG availability differ. Empty maps and PAT reads with no mapped CpG need review. | Boundary, chromosome-start, empty/no-overlap, duplicate interval, reverse-strand, insertion, and soft-clip cases with explicit coordinate expectations. |
| High: atlas inputs | Coordinate mismatches are printed and skipped; theta list parsing uses `eval()`. | Define expected mismatch handling; validate safe list parsing and schema failures with small synthetic tables. |
| Medium: error/output behavior | Missing CpG configuration calls `sys.exit()` in a module without importing `sys`; no output rows leave an existing output untouched. | Missing-input and no-coverage cases with explicit exceptions and output lifecycle expectations. |
| Medium: less-used helpers | SNP reader return shape differs from the shared reader contract; coordinate rows with no CpGs omit `origin`; `jsonconverter()` references unimported `datetime`. | Small direct API tests before relying on these paths. |
| Medium: packaging/integration | Sample JSON is invalid, example output paths target `/home/appuser`, and Docker expects external code. | Repair examples within a scoped task; make dependencies/build context explicit and verify a small CLI/container run. |
| Separate: upstream integration | Eight upstream commits followed the selected base; parser behavior overlaps local edits. | Review `git log 687b9d3..origin/master` and diffs; agree on merge/cherry-pick strategy and establish tests before integrating. Preserve the restoration commit. |

## Safe resume and maintenance

1. Read the project context map, then this plan and the architecture document for
   the intended change. Inspect `git status --short --branch` and current diffs.
2. Expect branch `ben-bam-support` and the original repository as `origin`.
   Verify current refs rather than assuming the dated snapshot above is still current.
3. Preserve any uncommitted documentation or user edits. Do not replay the Git
   restoration commands, replace local files from `master`, or remove backup data.
4. For runtime work, select an isolated environment and install the package there
   under that task's scope. A proposed baseline command is
   `PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider` from the repo
   root. It has not been executed successfully here. Check scientific outputs on
   representative fixtures before larger runs.
5. Update this plan with actual validation and the next safe step. Update
   architecture when behavior/contracts change and both indexes when active
   context files change. Keep unrelated repairs and upstream integration separate.
