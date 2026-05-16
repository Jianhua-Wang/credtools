# API Reference

The command line is the main interface, but CREDTOOLS also exposes Python
functions and objects. Use the API when you want to build notebooks, custom
workflows, or tests around the same logic used by the CLI.

Start with:

| Page | Contents |
| --- | --- |
| [Core Objects](core.md) | `Locus`, `LocusSet`, `LDMatrix`, summary-stat loading |
| [Preprocessing](preprocessing.md) | munging and chunking helpers |
| [Meta-Analysis](meta.md) | meta-analysis functions |
| [Quality Control](qc.md) | QC metrics and outlier removal |
| [Fine-Mapping](finemapping.md) | pipeline, fine-mapping, credible sets |
| [Plotting](plotting.md) | plotting functions |

!!! note "API stability"
    The CLI is the stable user surface. Internal helper functions may change
    faster than command behavior.
