# Interactive CLI plan

This document records the intended interactive CLI direction. It is not implemented yet.

## Decision

The target user-facing command surface is documented in `docs/CLI_INTERFACE.md`:

```bash
evolvepro init [DIR]
evolvepro evolve
evolvepro validate
evolvepro plot
```

`evolvepro init` owns interactive workspace setup. A separate public `configure` command is not part of the target UX; reconfiguration should go through `evolvepro init --reconfigure`.

Keep the current `argparse` CLI as the compatibility contract until the target commands have regression coverage:

```bash
evolvepro init [DIR]
evolvepro run --config evolvepro.toml
evolvepro validate [ROUND_DIR]
evolvepro process ...
evolvepro dms grid-search ...
```

Do not switch to `rich-click` until the workspace/config/lock/manifest contract is stable and public workflows have CLI regression coverage.

The future interactive CLI should use `rich-click` or an equivalent Click/Rich stack.

## Target behavior

Interactive behavior belongs primarily to `evolvepro init` and `evolvepro init --reconfigure`.

Interactive flow:

1. Detect whether the target workspace already exists.
2. Refuse to overwrite existing configuration unless `--reconfigure` is present.
3. If overriding existing configuration, prompt `Override it? [Y/n]` unless `--yes` is present.
4. Collect protein/workspace settings.
5. Write deterministic workspace config.
6. Initialize lock/index state.

Execution remains non-interactive. `evolvepro evolve` reads workspace state, validates previous artifacts, and advances the state machine.

## Source of truth

The source of truth remains file-based:

- `evolvepro.toml`
- `evolvepro.lock`
- `rounds/<round-id>/round.lock`
- `rounds/<round-id>/manifest.json`
- CAS snapshots under `.evolvepro/cas/`

Interactive prompts must write the same workspace config consumed by `evolvepro evolve`.

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
evolvepro init . --set-default --yes
evolvepro init . --reconfigure
evolvepro evolve
evolvepro validate
```

Interactive setup should call the same Python functions as non-interactive setup flags.

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
