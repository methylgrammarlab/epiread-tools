# epiread-tools project instructions

This Python package reads methylation epireads, PAT records, and modified-base
BAM files into sparse matrices, and produces per-CpG methylation summaries. It
also provides atlas readers used by external deconvolution software.

## Start here

Read [prompts/README.md](prompts/README.md) for the project context map.
Before implementation, validation, or maintenance, read
[prompts/plans.md](prompts/plans.md). Before changing architecture, interfaces,
formats, defaults, or important behavior, also read
[prompts/architecture.md](prompts/architecture.md).

## Active context index

| File | Purpose |
| --- | --- |
| [prompts/README.md](prompts/README.md) | Context routing, update rules, and archive index. |
| [prompts/architecture.md](prompts/architecture.md) | Observed architecture, data contracts, defaults, and limitations. |
| [prompts/plans.md](prompts/plans.md) | Current state, Git restoration record, validation, and follow-up work. |

Adding, renaming, or deleting an active prompt file requires updating this index
and `prompts/README.md`. Archive-only changes update the nearest archive index
unless an active route or dependency also changes.

## Working rules

- Inspect Git status before edits and preserve unrelated user changes. Keep
  changes focused; record a short plan and meaningful verification for substantial work.
- Keep architecture and planning documents current as implementation changes.
  Label proposals, uncertainties, and unrun checks explicitly.
- Human-authored `prompts/prompt_dev_*.md` files are read-only unless editing
  them is explicitly authorized. None exist at project setup.
- The active branch is `ben-bam-support`, deliberately based on upstream
  `687b9d3`. Read the restoration record before attempting upstream integration;
  do not replace this historical base as part of unrelated maintenance.
- Preserve the small fixtures and their indexes in `tests/data/`, including
  the currently unused variants. `.gitignore` excludes new generated artifacts;
  existing tracked `dist/` packages and `.idea/` settings remain tracked.
- For genomic changes, establish assembly and coordinate compatibility from
  input evidence. Preserve sorted/indexed inputs where readers require them;
  validate counts, boundaries, strand behavior, and output values on small fixtures.
  Existing mixed coordinate behavior is documented, not a uniform convention.
- Use the existing Python toolchain and tests. Do not infer correctness from
  syntax checks or the existing coverage-only BAM test. Record missing dependencies
  and functional checks that were not run.
- Commit, push, or publish within the user's authorized scope. Report unexpected
  Git obstacles before choosing a materially different integration strategy.

Global agent startup and private memory configuration live outside this scaffold.
