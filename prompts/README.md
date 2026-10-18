# Project context map

This directory holds the active technical context for epiread-tools. Start with
the project [AGENTS.md](../AGENTS.md), then read the documents relevant to the task.

| File | Read when | Contents |
| --- | --- | --- |
| [README.md](README.md) | Starting or updating project context. | Active file map, update rules, and deep-record index. |
| [architecture.md](architecture.md) | Learning the codebase or changing architecture, interfaces, formats, defaults, or behavior. | Modules, data flow, contracts, deployment, and observed limitations. |
| [plans.md](plans.md) | Implementing, validating, maintaining, or resuming work. | Current branch, restoration history, validation evidence, open issues, and next steps. |

The root [README.md](../README.md) is the existing user-facing guide;
[README.txt](../README.txt) contains historical operational commands with
machine-specific paths. Their configuration inconsistencies are recorded in the
active context; neither supersedes the observed implementation.

## Update rules

- Update `architecture.md` when the implemented design, interfaces, formats,
  defaults, dependencies, or important behavior change.
- Update `plans.md` after meaningful work with current state, validation results,
  remaining risks, and safe-resume instructions. Distinguish completed work from
  proposed follow-ups; a backlog entry is not authorization to implement it.
- Adding, renaming, or deleting active prompt files requires updating both this
  map and the index in `AGENTS.md`.
- Optional `prompt_dev_*.md` files are human-authored implementation prompts.
  Treat them as read-only unless editing is explicitly authorized. None currently exist.
- Dated `*_deep_records_*.md` files are historical context and are excluded from
  routine development reads. If active material is later archived, preserve its
  provenance and add a dated index entry describing source files and purpose.
- Archive-only changes update the nearest archive index; update active routing
  too if an active route or dependency changes.

## Deep Records Index

Deep records preserve historical material. Read one only when the user asks for
history or an active context file points to it for a specific current question.

- None yet.
