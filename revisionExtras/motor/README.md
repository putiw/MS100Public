# Reviewer-requested Motor analysis

This package is additive and self-contained. It uses the attached 13-tract
definition verbatim, records missing components, and never encodes a missing
component as a zero-valued metric. Subject-level masks and maps belong only in
the ignored local scratch directory configured by the user.

The master runner is `code/run_all_motor.py`. It requires a local staged input
manifest and aggregate metric tables; it refuses to publish model values when
those inputs are absent. It produces the complete 196-cell model grid with an
explicit status for every cell, so an unavailable input cannot be mistaken for
a silently omitted result.

No paired ΔAUC, DeLong test, ΔR², or new multiplicity family is implemented.

