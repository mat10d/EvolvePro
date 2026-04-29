# Interactive CLI plan

This document records the intended interactive CLI direction. It is not implemented yet.

## Decision

Keep the current `argparse` CLI as the non-interactive compatibility contract for now:

```bash
evolvepro init [DIR]
evolvepro run --config evolvepro.toml
evolvepro validate [ROUND_DIR]
```

Do not switch to `rich-click` until the workspace/config/lock/manifest contract is stable and public workflows have CLI regression coverage.

The future interactive CLI should use `rich-click` or an equivalent Click/Rich stack.

## Target behavior

Running `evolvepro` with no subcommand should eventually enter interactive mode.

Interactive flow:

1. Detect the current workspace or offer to initialize one.
2. Show a command menu.
3. Let the user choose a workflow.
4. Collect or edit workflow settings.
5. Write deterministic TOML.
6. Confirm before execution.
7. Execute through the same non-interactive run path.
8. Write lock and manifest files.

Interactive mode is a config generator/editor. It must never become a separate execution path.

## Source of truth

The source of truth remains file-based:

- `evolvepro.toml`
- `rounds/<round-id>/round.toml` when round-specific configs are introduced
- `rounds/<round-id>/round.lock.json`
- `rounds/<round-id>/manifest.json`

Interactive prompts must write the same TOML schema used by `evolvepro run --config`.

## Idempotency requirements

Interactive writes must be deterministic:

- stable TOML key order
- no timestamps
- no random UUIDs
- paths written relative to the config/workspace when possible
- same selections produce byte-identical TOML
- existing user config is not overwritten without explicit confirmation

## UX goals

Future interactive UX should provide:

- command selection
- workspace status display
- available round selection
- input file picker constrained to workspace paths where possible
- config diff preview before write
- validation status before run
- clear hash mismatch explanations
- explicit accept/update flow when input files changed

## Non-interactive parity

Every interactive action must map to a non-interactive command.

Examples:

```bash
evolvepro init .
evolvepro run --config evolvepro.toml
evolvepro validate outputs/kelsic/Round1
```

Interactive mode should call the same Python functions as these commands.

## Dependency policy

Do not add `rich-click` yet.

Add it only when implementing interactive scaffolding. When added, it becomes a runtime dependency and must be locked in `uv.lock`.

## Migration order

1. Keep current `argparse` CLI stable for non-interactive tests.
2. Finish README/script workflow migration to CLI commands.
3. Stabilize config and round layout based on those migrations.
4. Add `rich-click` dependency.
5. Replace or wrap the parser with Click/Rich command definitions.
6. Add interactive-only tests for prompt routing and TOML generation.
7. Keep non-interactive CLI tests passing unchanged.
