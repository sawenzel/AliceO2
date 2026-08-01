"""Stream A — the B-rep -> CSG simplification pipeline.

See `scripts/geometry/CSG_Pipeline.md` for the brief and `scripts/geometry/Stream_A_CSG.md`
for the measurements this package produced.

Modules
-------
`occ_env`   re-exec into the alibuild Python 3.10 that can import pythonOCC.
`census`    step 1 of the plan: measure, per solid, what representation could possibly apply.
"""
