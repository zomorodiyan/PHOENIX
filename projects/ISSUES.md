# Issues

Internal bug and correctness tracker. Each issue lives in its own folder alongside project folders, following the same naming convention:

- Open:    `YYYYMMDD_ISSUE-XXX_TITLE/`
- Fixed:   `YYYYMMDD_ISSUE-XXX_TITLE-FIXEDYYYYMMDD/`

See [ISSUE_TEMPLATE.md](ISSUE_TEMPLATE.md) for the issue format.

---

| ID        | Title                                                         | Severity | Status   | Discovered |
|-----------|---------------------------------------------------------------|----------|----------|------------|
| [ISSUE-001](20260411_ISSUE-001_DIRECTIONAL_SYMMETRY-FIXED20260413/) | `ishift_field` incorrect when `dj > 0`                        | Medium   | Fixed    | 2026-04-11 |
| [ISSUE-002](20260413_ISSUE-002_TOP_SURFACE_HEAT_LOSS/)              | Top surface convective loss missing; radiation double-counted | High     | Open     | 2026-04-13 |
| [ISSUE-003](20260413_ISSUE-003_BOTTOM_SURFACE_HLOSSCONVEC/)         | Bottom surface `hlossconvec` silently includes radiation      | Medium   | Open     | 2026-04-13 |
