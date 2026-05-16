# CLI Reference

This section is for lookup. If you are learning the workflow, use the tutorials
first.

## Commands

| Command | Use it when |
| --- | --- |
| `credtools munge` | raw summary statistics need cleaning and standard column names |
| `credtools chunk` | genome-wide files need to be split into loci |
| `credtools meta` | you want to run meta-analysis as a separate step |
| `credtools qc` | you want standalone QC tables or cleaned inputs |
| `credtools finemap` | you want standalone fine-mapping |
| `credtools pipeline` | you want meta-analysis, QC, and fine-mapping in one run |
| `credtools plot` | you want quick QC or result plots |

General flags:

```bash
credtools --help
credtools --version
credtools --log-file credtools.log
```

Command-specific help:

```bash
credtools pipeline --help
credtools finemap --help
```

!!! tip "Use the long flag in scripts"
    Short flags are handy at the terminal. Long flags make scripts easier to
    read six months later.
